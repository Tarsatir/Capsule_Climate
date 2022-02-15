using Printf
using Statistics
using Distributions
using StatsBase
using Random
using Agents
using BenchmarkTools
using TimerOutputs
using RecursiveArrayTools

# Include files
include("../results/write_results.jl")
include("custom_schedulers.jl")
include("global_parameters.jl")
include("misc.jl")
include("government.jl")
include("labormarket.jl")
include("consumermarket.jl")
include("household.jl")
include("consumer_good.jl")
include("capital_good.jl")
include("general_producers.jl")
include("macro.jl")


"""
Initializes model.
- Receives:
    * number of agents of each type 
    * number of initial employees per producer type.
- Returns:
    * model: Agent Based Model struct
    * global_param: struct containing global parameter values
    * macro_struct: mutable struct containing macro variables
    * gov_struct: mutable struct containing variables concerning the government
    * labormarket_struct: mutable struct containing variables concerning the labor market
"""
function initialize_model(;
    n_captlgood = 50,
    n_consrgood = 200,
    n_households = 2500,
    n_init_emp_cp = 11,
    n_init_emp_kp = 3,
    n_bp = 7,
    n_lp = 7
    )

    # Initialise model struct
    model = ABM(Union{Household, CapitalGoodProducer, ConsumerGoodProducer};
                scheduler = by_type(true,true))

    # Initialize struct that holds global params
    global_param  = initialize_global_params()

    # Initialize struct that holds macro variables
    macro_struct = initialize_macro()

    # Initialize labor market struct
    labormarket_struct = initialize_labormarket()

    # Initialize government struct
    gov_struct = initialize_government()

    # Global id
    id = 1

    # Initialize households
    for hh_i in 1:n_households

        hh = initialize_hh(id)
        add_agent!(hh, model)

        id += 1
    end

    # Initialize consumer good producers
    for cp_i in 1:n_consrgood

        # Decide if producer makes basic or luxury goods
        # In init, half of producers are allocated basic and half luxury
        type_good = "Basic"
        if cp_i > n_consrgood / 2
            type_good = "Luxury"
        end

        # Initialize capital good stock
        # Machines have random age as to allow replacement in early periods
        machine_struct = initialize_machine(global_param.η)

        cp = initialize_cp(id, machine_struct, n_consrgood, type_good, n_init_emp_cp)
        add_agent!(cp, model)

        id += 1
    end

    # Initialize capital good producers
    for kp_i in 1:n_captlgood

        kp = initialize_kp(id, kp_i, n_captlgood, n_init_emp_kp)
        add_agent!(kp, model)

        id += 1
    end

    # Initialize schedulers
    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = schedule_per_type(true, model)

    # Let all capital good producers select historical clients
    for kp_id in all_kp
        select_HC_kp!(model[kp_id], all_cp)
    end

    # Let all households select cp of both types for trading network
    for hh_id in all_hh
        select_bp_lp_hh!(model[hh_id], all_bp, all_lp, n_bp, n_lp)
    end

    # Spread employed households over producerss
    spread_employees_lm!(
        labormarket_struct,
        all_hh, 
        all_cp, 
        all_kp,
        n_init_emp_cp,
        n_init_emp_kp,
        model
    )

    # Update market shares of cp
    update_marketshares_cm!(all_cp, model)

    return model, global_param, macro_struct, gov_struct, labormarket_struct
end


function model_step!(
    global_param, 
    macro_struct, 
    gov_struct, 
    labormarket_struct,
    model::ABM
    )

    # Update schedulers
    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = schedule_per_type(true, model)

    # Reset kp brochures of all cp
    for cp_id in all_cp
        reset_brochures_cp!(model[cp_id])
    end

    # (1) capital good producers innovate and send brochures

    # Determine distance matrix between capital good producers
    kp_distance_matrix = get_capgood_euclidian(all_kp, model)

    for kp_id in all_kp
        kp = model[kp_id]
        reset_order_queue_kp!(kp)
        innovate_kp!(kp, global_param, all_kp, kp_distance_matrix, model)
        send_brochures_kp!(kp, all_cp, global_param, model)
    end

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    for cp_id in all_cp
        cp = model[cp_id]
        plan_production_cp!(cp, global_param, model)
        plan_investment_cp!(cp, global_param, all_kp, model)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp_id in all_kp
        plan_production_kp!(model[kp_id])
    end

    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))
    # println(length(labormarket_struct.employed))

    # (3) labor market matching process
    labormarket_process!(
        labormarket_struct,
        all_hh, 
        all_p,
        global_param.ϵ,
        gov_struct.UB,
        model
    )
    update_avg_T_unemp_lm(labormarket_struct, model)


    # (4) Producers pay workers their wage. Government pays unemployment benefits
    for p_id in all_p
        pay_workers_p!(model[p_id], model)
    end
    pay_unemployment_benefits_gov!(gov_struct, labormarket_struct.unemployed, model)

    # Government receives income taxes
    levy_income_tax_gov!(gov_struct, all_hh, model)


    # (5) Production takes place for cp and kp
    for cp_id in all_cp
        produce_goods_cp!(model[cp_id])
    end

    for kp_id in all_kp
        produce_goods_kp!(model[kp_id])
    end

    
    # (6) Transactions take place on consumer market

    # Households update wealth level
    for hh_id in all_hh
        update_wealth_hh!(model[hh_id])
    end

    # Consumer market process
    consumermarket_process!(
        all_hh,
        all_cp,
        all_bp,
        all_lp,
        gov_struct,
        model
    )

    # cp make up profits
    for cp_id in all_cp
        compute_profits_cp!(model[cp_id])
    end

    # (6) kp deliver goods to cp, kp make up profits
    for kp_id in all_kp
        send_orders_kp!(model[kp_id], model)
    end

    for cp_id in all_cp
        increase_machine_age_cp!(model[cp_id])
    end

    # (7) government receives profit taxes and computes budget balance
    levy_profit_tax_gov!(gov_struct, all_p, model)
    compute_budget_balance(gov_struct)

    # (7) macro-economic indicators are updated.
    update_macro_timeseries(
        macro_struct, 
        all_hh, 
        all_cp, 
        all_kp,
        labormarket_struct.E,
        gov_struct.curr_acc.Exp_UB[end],
        model
    )

    # TODO update market share cp

    # Remove bankrupt companies and introduce new companies.

    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
end

to = TimerOutput()

# df_agent = init_agent_dataframe(model)

@timeit to "init" model, global_param, macro_struct, gov_struct, labormarket_struct = initialize_model()
for i in 1:100
    println("Step ", i)
    @timeit to "step" model_step!(global_param, macro_struct, gov_struct, labormarket_struct, model)
end

# @timeit to "step" run!(model, dummystep, model_step!, 10)

println(macro_struct.GDP)

@timeit to "save macro" save_macro_data(macro_struct)

all_hh, all_cp, all_kp, all_bp, all_lp, all_p = schedule_per_type(true, model)

@timeit to "save findist" save_final_dist(all_hh, model)

show(to)
println()