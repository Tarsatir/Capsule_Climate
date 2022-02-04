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
include("model.jl")
include("misc.jl")
include("government.jl")
include("labormarket.jl")
include("consumermarket.jl")
include("household.jl")
include("consumer_good.jl")
include("capital_good.jl")
include("general_producers.jl")
include("macro.jl")



function initialize_model(;
    n_captlgood = 50,
    n_consrgood = 200,
    n_households = 2500,
    n_init_emp_cp = 11,
    n_init_emp_kp = 3
    )

    # initialise model struct
    model = ABM(Union{Household, CapitalGoodProducer, ConsumerGoodProducer};
                scheduler = by_type(true,true))

    # initialize struct that holds global params
    global_param  = initialize_global_params()

    # initialize struct that holds macro variables
    macro_struct = initialize_macro()

    # initialize labor market struct
    labormarket_struct = initialize_labormarket()

    # initialize consumer market struct
    # consumermarket_struct = initialize_consumermarket()

    # initialize government struct
    gov_struct = initialize_government()

    # global id
    id = 1

    # initialize households
    for hh_i in 1:n_households

        hh = initialize_hh(id)
        add_agent!(hh, model)

        id += 1
    end

    # initialize consumer good producers
    for cp_i in 1:n_consrgood

        # decide if producer makes basic or luxury goods
        type_good = "Basic"
        if cp_i > n_consrgood / 2
            type_good = "Luxury"
        end

        # initialize capital good stock
        machine_struct = initialize_machine()
        machine_struct.age = rand(0:global_param.η)

        cp = initialize_cp(id, machine_struct, n_consrgood, type_good, n_init_emp_cp)

        # push!(all_agents.all_cp, cp)
        
        # if type_good == "Basic"
        #     push!(all_agents.all_bp, cp.cp_id)
        # else
        #     push!(all_agents.all_lp, cp.cp_id)
        # end

        add_agent!(cp, model)
        id += 1
    end

    # initialize capital good producers
    for kp_i in 1:n_captlgood

        # make choice for historical clients
        # HC = sample(all_cp, 10; replace=false)

        # initialize capital good producer
        kp = initialize_kp(id, kp_i, n_captlgood, n_init_emp_kp)

        # push!(all_agents.all_kp, kp)
        add_agent!(kp, model)
        id += 1
    end

    # initialize schedulers
    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = per_type(true, model)

    for kp_id in all_kp
        kp = model[kp_id]
        select_HC_kp!(kp, all_cp)
    end

    # spread employed households over producerss
    spread_employees_lm!(
        labormarket_struct,
        all_hh, 
        all_cp, 
        all_kp,
        n_init_emp_cp,
        n_init_emp_kp,
        model
    )

    # return model, all_agents, global_param, macro_struct, gov_struct, labormarket_struct, consumermarket_struct
    return model, global_param, macro_struct, gov_struct, labormarket_struct
end


function model_step!(model, 
    global_param, 
    macro_struct, 
    gov_struct, 
    labormarket_struct
    )

    # update schedulers
    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = per_type(true, model)

    # reset brochures of all consumer good producers
    for cp_id in all_cp
        cp = model[cp_id]
        reset_brochures_cp!(cp)
    end

    # (1) capital good producers innovate and send brochures

    # determine distance matrix between capital good producers
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
        plan_production_cp!(cp, global_param)
        plan_investment_cp!(cp, global_param, all_kp, model)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp_id in all_kp
        kp = model[kp_id]
        plan_production_kp!(kp)
    end

    println(sum(map(p_id -> length(model[p_id].Emp), all_p)))
    println(length(labormarket_struct.employed))

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


    # (5) Production takes place for cp and kp
    for cp_id in all_cp
        produce_goods_cp!(model[cp_id])
    end

    for kp_id in all_kp
        produce_goods_kp!(model[kp_id])
    end

    # Government receives income taxes
    levy_income_tax_gov!(gov_struct, all_hh, model)
    # compute_budget_balance(gov_struct)

    # (6) Households pick prefered products to buy and set budget and consumption package
    for hh_id in all_hh
        compute_exp_income_hh!(model[hh_id], 
                               labormarket_struct.P_HU, 
                               labormarket_struct.P_UU, 
                               gov_struct.UB,
                               model)
        set_savingsrate_hh!(model[hh_id], labormarket_struct.avg_T_unemp, gov_struct.UB)
    end


    # (6) Transactions take place on consumer market
    consumermarket_process!(
                            # consumermarket_struct,
                            all_hh,
                            all_cp,
                            all_bp,
                            all_lp,
                            gov_struct,
                            model)

    # cp make up profits
    for cp_id in all_cp
        compute_profits_cp!(model[cp_id])
    end

    # (6) kp deliver goods to cp, kp make up profits
    for kp_id in all_kp
        send_orders_kp!(model[kp_id], model)
    end

    # (7) government receives profit taxes
    # TODO

    # (7) macro-economic indicators are updated.
    update_macro_stats(macro_struct, 
                       all_hh, 
                       all_cp, 
                       all_kp,
                       labormarket_struct.E,
                       gov_struct.curr_acc.Exp_UB[end],
                       model)

    # TODO update market share cp

    println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))


end

to = TimerOutput()

@timeit to "init" model, global_param, macro_struct, gov_struct, labormarket_struct = initialize_model()
for i in 1:100
    println("Step ", i)
    @timeit to "step" model_step!(model, global_param, macro_struct, gov_struct, labormarket_struct)
end

println(macro_struct.GDP)

@timeit to "save" save_macro_data(macro_struct)

show(to)
println()