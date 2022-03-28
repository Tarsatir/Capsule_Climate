using Printf
using Statistics
using Distributions
using StatsBase
using Random
using Agents
using BenchmarkTools
using TimerOutputs
using RecursiveArrayTools
# using PyCall

# Include files
include("../results/write_results.jl")
include("helpers/custom_schedulers.jl")
include("helpers/dist_matrix.jl")
include("global_parameters.jl")

include("objects/accounting_firms.jl")
include("objects/accounting_govt.jl")
include("objects/machine.jl")

include("agents/government.jl")
include("agents/household.jl")
include("agents/consumer_good.jl")
include("agents/capital_good.jl")
include("agents/general_producers.jl")
include("agents/bank.jl")
include("agents/indexfund.jl")

include("macro_markets/labormarket.jl")
include("macro_markets/consumermarket.jl")
include("macro_markets/macro.jl")

# TEMP
include("../plotting/plot_D_p.jl")


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
function initialize_model(
    T;
    n_captlgood = 50,
    n_consrgood = 200,
    n_households = 2500,
    n_init_emp_cp = 10,
    n_init_emp_kp = 4,
    n_bp = 7,
    n_lp = 7,
    changed_params
    )

    # Initialise model struct
    model = ABM(Union{Household, CapitalGoodProducer, ConsumerGoodProducer};
                scheduler = by_type(true,true),
                warn=false)

    # Initialize struct that holds global params
    global_param  = initialize_global_params(changed_params)

    # Initialize struct that holds macro variables
    macro_struct = initialize_macro(T)

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
        n_machines_init = 40  #TODO: put this in parameters

        machines = Vector{Machine}()
        K = n_init_emp_cp * 100
        for _ in 1:n_machines_init
            # Machines have random age as to allow replacement in early periods
            freq = K/n_machines_init
            machine_struct = initialize_machine(freq; η=global_param.η)
            push!(machines, machine_struct)
        end

        cp = initialize_cp(
                id, 
                machines,  
                type_good, 
                n_init_emp_cp, 
                global_param.μ1,
                global_param.ι;
                n_consrgood=n_consrgood
            )
        update_n_machines_cp!(cp)
        add_agent!(cp, model)

        id += 1
    end

    # Initialize capital good producers
    for kp_i in 1:n_captlgood

        kp = initialize_kp(id, kp_i, n_captlgood)
        add_agent!(kp, model)

        id += 1
    end

    # Initialize schedulers
    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = schedule_per_type(true, model)

    # Initialize bank struct
    bank_struct = initialize_bank(all_p, model)

    # Initialize index fund struct
    indexfund_struct = initialize_indexfund()

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

    for p_id in all_p
        close_balance_p!(
            model[p_id], 
            global_param.Λ,
            global_param.r,
            global_param.η,
            gov_struct.τᴾ,
            indexfund_struct
        )
    end

    return model, global_param, macro_struct, gov_struct, labormarket_struct, bank_struct, indexfund_struct
end


function model_step!(
    t::Int,
    T::Int,
    global_param::GlobalParam, 
    macro_struct::MacroEconomy, 
    gov_struct::Government, 
    labormarket_struct::LaborMarket,
    bank_struct::Bank,
    indexfund_struct::IndexFund,
    model::ABM
    )

    # Update schedulers
    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = schedule_per_type(true, model)

    # Update firm age
    for p_id in all_p
        model[p_id].age += 1
    end

    # Check if households still have enough producers, otherwise sample more
    refill_suppliers_all_hh!(
        all_hh,
        all_bp, 
        all_lp,
        model
    )

    # Clear current account, decide how many debts to repay, reset kp brochures of all cp
    for cp_id in all_cp
        clear_firm_currentaccount_p!(model[cp_id].curracc)
        reset_brochures_cp!(model[cp_id])
    end

    # (1) capital good producers innovate and send brochures

    # Determine distance matrix between capital good producers
    kp_distance_matrix = get_capgood_euclidian(all_kp, model)

    for kp_id in all_kp

        kp = model[kp_id]
        clear_firm_currentaccount_p!(model[kp_id].curracc)

        if length(macro_struct.w̄_avg) == 0
            w̄ = 1.0
        else
            w̄ = macro_struct.w̄_avg[end]
        end

        innovate_kp!(
            kp, 
            global_param, 
            all_kp, 
            kp_distance_matrix,
            w̄,
            model
        )
        send_brochures_kp!(kp, all_cp, global_param, model)
    end

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    for cp_id in all_cp
        cp = model[cp_id]
        
        # Plan production for this period
        plan_production_cp!(cp, global_param, model)

        # Plan investments for this period
        plan_investment_cp!(cp, global_param, all_kp, model)

        # See if enough funds available for investments and production, otherwise
        # change investments and production to match funding availability.
        check_funding_restrictions_cp!(cp, global_param.Λ, global_param.r, global_param.ωW, model)

        # Send orders to kp
        order_machines_cp!(cp, model)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp_id in all_kp
        plan_production_kp!(model[kp_id], model)
    end

    # (3) labor market matching process
    labormarket_process!(
        labormarket_struct,
        all_hh, 
        all_p,
        global_param.ϵ,
        global_param.max_g_wᴼ,
        gov_struct.UB,
        global_param,
        model
    )
    # TODO: check if this is still needed, otherwise delete
    # update_avg_T_unemp_lm(labormarket_struct, model)


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
    distribute_dividends_if!(
        indexfund_struct,
        all_hh,
        model
    )

    for hh_id in all_hh
        # Reset unsatisfied demand
        model[hh_id].unsat_dem = Vector()

        # Update household wealth
        update_wealth_hh!(model[hh_id])
    end

    # TODO: put this somewhere else
    # all_W_hh = map(hh_id -> model[hh_id].W[end], all_hh)
    all_W_hh = Vector{Float64}()

    # Consumer market process
    consumermarket_process!(
        all_hh,
        all_cp,
        all_bp,
        all_lp,
        all_W_hh,
        gov_struct,
        global_param,
        model
    )

    # Households decide to switch suppliers based on satisfied demand and prices
    for hh_id in all_hh
        decide_switching_hh!(
            model[hh_id],
            global_param.ψ_Q,
            global_param.ψ_P,
            all_bp,
            all_lp,
            model
        )
    end

    # (6) kp deliver goods to cp, kp make up profits
    for kp_id in all_kp
        send_orders_kp!(model[kp_id], model)
    end

    # Close balances of all firms
    for p_id in all_p
        # Update amount of owned capital, increase machine age
        if typeof(model[p_id]) == ConsumerGoodProducer
            update_n_machines_cp!(model[p_id])
            increase_machine_age_cp!(model[p_id])
        end

        # Close balances of firms, if insolvent, liquidate firms
        close_balance_p!(
            model[p_id], 
            global_param.Λ,
            global_param.r,
            global_param.η,
            gov_struct.τᴾ,
            indexfund_struct
        )
    end 

    # (7) government receives profit taxes and computes budget balance
    levy_profit_tax_gov!(gov_struct, all_p, model)
    compute_budget_balance(gov_struct)
    redistribute_surplus_gov!(gov_struct, all_hh, model)

    # Update market shares of cp
    update_marketshare_cp!(all_bp, all_lp, model)
    update_marketshare_kp!(all_kp, model)

    # Select producers that will be declared bankrupt and removed
    bankrupt_bp, bankrupt_lp, bankrupt_kp, bankrupt_kp_i = check_bankrupty_all_p!(all_p, all_kp, model)

    println("avg age bankrupt kp: $(mean(map(kp_id -> model[kp_id].age, bankrupt_kp)))")

    # (7) macro-economic indicators are updated.
    update_macro_timeseries(
        macro_struct, 
        all_hh, 
        all_cp, 
        all_kp,
        all_bp,
        all_lp,
        bankrupt_bp,
        bankrupt_lp,
        bankrupt_kp,
        labormarket_struct.E,
        gov_struct,
        global_param,
        model
    )

    # Remove bankrupt companies.
    kill_all_bankrupt_p!(
        bankrupt_bp, 
        bankrupt_lp, 
        bankrupt_kp, 
        all_hh,
        all_kp,
        labormarket_struct,
        indexfund_struct, 
        model
    )
    update_unemploymentrate_lm!(labormarket_struct)
    # println("E 4: ", labormarket_struct.E[end])

    # println("Number of cp: ", length(all_cp), ", Number of kp: ", length(all_kp))


    # Replace bankrupt kp. Do this before you replace cp, such that new cp can also choose
    # from new kp 
    replace_bankrupt_kp!(
        bankrupt_kp, 
        bankrupt_kp_i, 
        all_kp,
        global_param.φ3,
        global_param.φ4, 
        global_param.α2,
        global_param.β2,
        model
    )

    # Replace bankrupt companies with new companies
    replace_bankrupt_cp!(
        bankrupt_bp, 
        bankrupt_lp,
        bankrupt_kp, 
        all_hh,
        all_bp,
        all_lp,
        all_kp,
        macro_struct.cu[end],
        global_param.φ1,
        global_param.φ2,
        global_param.φ3,
        global_param.φ4,
        macro_struct.p̄[end],
        macro_struct.w̄_avg[end],
        global_param.μ1,
        global_param.ι, 
        model
    )

    # if t == T
    #     plot_D_p(all_bp, model)
    # end

    # println("All employed: ", sum(map(p_id -> length(model[p_id].employees), all_p)))

end


function run_simulation(;
    T=400::Int,
    changed_params=nothing,
    full_output=true::Bool
    )::Float64

    to = TimerOutput()

    @timeit to "init" model, global_param, macro_struct, gov_struct, labormarket_struct, bank_struct, indexfund_struct = initialize_model(T; changed_params=changed_params)
    for t in 1:T
        
        if t % 100 == 0
            println("Step $t")
        end

        @timeit to "step" model_step!(
                                t,
                                T, 
                                global_param, 
                                macro_struct, 
                                gov_struct, 
                                labormarket_struct, 
                                bank_struct, 
                                indexfund_struct,
                                model
                            )
    end

    @timeit to "save macro" save_macro_data(macro_struct)

    all_hh, all_cp, all_kp, all_bp, all_lp, all_p = schedule_per_type(true, model)

    @timeit to "save findist" save_final_dist(all_hh, all_bp, all_lp, all_kp, model)

    if full_output
        show(to)
        println()
    end

    return macro_struct.GDP[end]
end

run_simulation()