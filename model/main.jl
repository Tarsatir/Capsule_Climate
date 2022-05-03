using Printf
using Statistics
using Distributions
using StatsBase
using Random
using Agents
using BenchmarkTools
using TimerOutputs
using RecursiveArrayTools
using DataStructures
using Parameters
using Dates

# Include files
include("../results/write_results.jl")
include("helpers/custom_schedulers.jl")
include("helpers/dist_matrix.jl")
include("helpers/intermediate_data_storage.jl")
include("helpers/update.jl")
include("global_parameters.jl")
include("init_parameters.jl")

include("objects/accounting_firms.jl")
include("objects/accounting_govt.jl")
include("objects/machine.jl")
include("objects/powerplant.jl")
include("objects/climate.jl")

include("agents/government.jl")
include("agents/indexfund.jl")
include("macro_markets/macro.jl")
include("agents/household.jl")
include("agents/consumer_good_producer.jl")
include("agents/capital_good_producer.jl")
include("agents/general_producer.jl")
include("agents/energy_producer.jl")
include("agents/bank.jl")

include("macro_markets/labormarket.jl")
include("macro_markets/consumermarket.jl")


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
    T::Int,
    labormarket_is_fordist::Bool;
    changed_params
    )

    # Initialise model struct
    model = ABM(Union{Household, CapitalGoodProducer, ConsumerGoodProducer};
                scheduler = by_type(true,true), warn=false)

    # Initialize struct that holds global params and initial parameters
    global_param  = initialize_global_params(labormarket_is_fordist, changed_params)
    init_param = InitParam()

    # Initialize struct that holds macro variables
    macro_struct = MacroEconomy(T=T)

    # Initialize labor market struct
    labormarket_struct = LaborMarket()

    # Initialize government struct
    gov_struct = Government(curracc = GovCurrentAccount(T=T))

    # Initialize energy producer
    ep = initialize_energy_producer(T, init_param, global_param)

    # Initialize climate struct
    climate_struct = Climate(T=T)

    # Global id
    id = 1

    # Initialize households
    hh_skills = sample_skills_hh(init_param)
    for i in 1:init_param.n_hh

        hh = Household(id=id, skill=hh_skills[i])
        hh.wʳ = max(gov_struct.w_min, hh.wʳ)
        add_agent!(hh, model)

        id += 1
    end

    # Initialize consumer good producers
    for _ in 1:init_param.n_cp

        # Initialize capital good stock
        machines = initialize_machine_stock(global_param.freq_per_machine, 
                        init_param.n_machines_init,
                        η=global_param.η; 
                        A_LP = init_param.A_LP_0,
                        A_EE = init_param.A_EE_0,
                        A_EF = init_param.A_EF_0
                   )

        cp = initialize_cp(
                id,
                machines,  
                init_param.n_init_emp_cp, 
                global_param.μ1,
                global_param.ι;
                n_consrgood=init_param.n_cp
            )
        update_n_machines_cp!(cp, global_param.freq_per_machine)
        add_agent!(cp, model)

        id += 1
    end

    # Initialize capital good producers
    for kp_i in 1:init_param.n_kp

        kp = initialize_kp(
                id, 
                kp_i, 
                init_param.n_kp; 
                A_LP=init_param.A_LP_0,
                A_EE=init_param.A_EE_0,
                A_EF=init_param.A_EF_0, 
                B_LP=init_param.B_LP_0,
                B_EE=init_param.B_EE_0,
                B_EF=init_param.B_EF_0
             )
        add_agent!(kp, model)

        id += 1
    end

    # Initialize schedulers
    all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

    # Initialize index fund struct
    indexfund_struct = IndexFund(T=T)

    # Let all capital good producers select historical clients
    for kp_id in all_kp
        select_HC_kp!(model[kp_id], all_cp)
    end

    # Let all households select cp of both types for trading network
    for hh_id in all_hh
        select_cp_hh!(model[hh_id], all_cp, init_param.n_cp_hh)
    end

    # Spread employed households over producerss
    spread_employees_lm!(
        labormarket_struct,
        gov_struct,
        all_hh, 
        all_cp, 
        all_kp,
        init_param.n_init_emp_cp,
        init_param.n_init_emp_kp,
        model
    )

    for p_id in all_p
        update_mean_skill_p!(model[p_id], model)
    end

    close_balance_all_p!(all_p, global_param, gov_struct.τᴾ,
                         indexfund_struct, 0, model)

    cm_dat = CMData(n_hh=init_param.n_hh, n_cp=init_param.n_cp)

    return model, global_param, init_param, macro_struct, gov_struct, ep, labormarket_struct, indexfund_struct, climate_struct, cm_dat
end


function model_step!(
    t::Int,
    to,
    global_param::GlobalParam, 
    init_param::InitParam,
    macro_struct::MacroEconomy, 
    gov_struct::Government,
    ep::EnergyProducer, 
    labormarket_struct::LaborMarket,
    indexfund_struct::IndexFund,
    climate_struct::Climate,
    cm_dat::CMData,
    model::ABM
    )

    # Update schedulers
    @timeit to "schedule" all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

    # Update firm age
    for p_id in all_p
        model[p_id].age += 1
    end

    # Check if households still have enough bp and lp, otherwise sample more
    refill_suppliers_all_hh!(
        all_hh,
        all_cp,
        init_param.n_cp_hh,
        model
    )

    # Clear current account, decide how many debts to repay, reset kp brochures of all cp
    for cp_id in all_cp
        model[cp_id].curracc = clear_firm_currentaccount_p!(model[cp_id].curracc)
        reset_brochures_cp!(model[cp_id])
    end

    # (1) kp and ep innovate, and kp send brochures

    # Determine distance matrix between kp
    @timeit to "dist mat" kp_distance_matrix = get_capgood_euclidian(all_kp, model)

    for kp_id in all_kp

        model[kp_id].curracc = clear_firm_currentaccount_p!(model[kp_id].curracc)

        @timeit to "kp innov" innovate_kp!(
            model[kp_id], 
            global_param, 
            all_kp, 
            kp_distance_matrix,
            macro_struct.w̄_avg[max(t-1,1)],
            t,
            ep,
            model
        )
        send_brochures_kp!(model[kp_id], all_cp, global_param, model)
    end

    # ep innovate
    innovate_ep!(ep, global_param, t)

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    for cp_id in all_cp
        cp = model[cp_id]
        
        # Plan production for this period
        @timeit to "plan prod cp" plan_production_cp!(cp, global_param, t, model)

        # Plan investments for this period
        @timeit to "plan inv cp" plan_investment_cp!(cp, all_kp, global_param, ep, t, model)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp_id in all_kp
        @timeit to "plan prod kp" plan_production_kp!(model[kp_id], global_param, model)
    end

    # (3) labor market matching process
    @timeit to "labormarket" labormarket_process!(
        labormarket_struct,
        all_hh, 
        all_p,
        global_param,
        gov_struct,
        model
    )

    # Update mean skill level of employees
    for p_id in all_p
        update_mean_skill_p!(model[p_id], model)
    end


    # (4) Producers pay workers their wage. Government pays unemployment benefits
    for p_id in all_p
        pay_workers_p!(model[p_id], model)
    end
    pay_unemployment_benefits_gov!(gov_struct, labormarket_struct.unemployed, t, model)

    # Government receives income taxes
    levy_income_tax_gov!(gov_struct, all_hh, t, model)


    # (5) Production takes place for cp and kp

    # Update energy prices
    compute_pₑ_ep!(ep, t)

    @timeit to "consumer prod" for cp_id in all_cp
        produce_goods_cp!(model[cp_id], ep, global_param, t)
    end

    @timeit to "capital prod" for kp_id in all_kp
        produce_goods_kp!(model[kp_id], ep, global_param, t)
    end

    # Let energy producer meet energy demand
    produce_energy_ep!(ep, all_cp, all_kp, global_param, indexfund_struct, t, model)

    
    # (6) Transactions take place on consumer market

    # Households update wealth level
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
    @timeit to "consumermarket" consumermarket_process!(
                                    all_hh,
                                    all_cp,
                                    all_W_hh,
                                    gov_struct,
                                    global_param,
                                    cm_dat,
                                    t,
                                    model,
                                    to
                                )

    # Households decide to switch suppliers based on satisfied demand and prices
    @timeit to "decide switching hh" decide_switching_all_hh!(
        global_param,
        all_hh,
        all_cp,
        all_p,
        init_param.n_cp_hh,
        model
    )

    # (6) kp deliver goods to cp, kp make up profits
    for kp_id in all_kp
        send_ordered_machines_kp!(model[kp_id], global_param, model)
    end

    # Close balances of all firms
    for cp_id in all_cp
        # Update amount of owned capital, increase machine age
        update_n_machines_cp!(model[cp_id], global_param.freq_per_machine)
        increase_machine_age_cp!(model[cp_id])
    end 

    # Close balances of firms, if insolvent, liquidate firms
    @timeit to "close balance" close_balance_all_p!(all_p, global_param, gov_struct.τᴾ,
                                    indexfund_struct, t, model)

    # (7) government receives profit taxes and computes budget balance
    levy_profit_tax_gov!(gov_struct, all_p, t, model)
    compute_budget_balance(gov_struct, t)
    redistribute_surplus_gov!(gov_struct, all_hh, model)

    # Update market shares of cp
    update_marketshare_p!(all_cp, model)
    update_marketshare_p!(all_kp, model)

    # Select producers that will be declared bankrupt and removed
    # @timeit to "check bankr" bankrupt_bp, bankrupt_lp, bankrupt_kp, bankrupt_kp_i = check_bankrupty_all_p!(all_p, all_kp, global_param, model)
    @timeit to "check bankr" bankrupt_cp, bankrupt_kp, bankrupt_kp_i = check_bankrupty_all_p!(all_p, all_kp, global_param, model)


    # (7) macro-economic indicators are updated.
    @timeit to "update macro ts" update_macro_timeseries(
        macro_struct,
        t, 
        all_hh, 
        all_cp, 
        all_kp,
        all_p,
        ep,
        bankrupt_cp,
        bankrupt_kp,
        labormarket_struct,
        gov_struct,
        indexfund_struct,
        global_param,
        model
    )

    # Update climate parameters, compute new carbon equilibria and temperature change
    collect_emissions_cl!(climate_struct, all_cp, all_kp, ep, t, model)
    carbon_equilibrium_tempchange_cl!(climate_struct, t)

    # Remove bankrupt companies.
    @timeit to "kill bankr p" kill_all_bankrupt_p!(
        bankrupt_cp,
        bankrupt_kp, 
        all_hh,
        all_kp,
        labormarket_struct,
        indexfund_struct, 
        model
    )
    update_unemploymentrate_lm!(labormarket_struct)


    # Replace bankrupt kp. Do this before you replace cp, such that new cp can also choose
    # from new kp
    @timeit to "replace kp" replace_bankrupt_kp!(
        bankrupt_kp, 
        bankrupt_kp_i, 
        all_kp,
        global_param,
        indexfund_struct,
        init_param,
        macro_struct,
        t,
        model
    )

    # Replace bankrupt companies with new companies
    @timeit to "replace cp" replace_bankrupt_cp!(
        bankrupt_cp, 
        bankrupt_kp, 
        all_hh,
        all_cp,
        all_kp,
        global_param,
        indexfund_struct,
        macro_struct,
        t,
        model
    )

    # Redistrubute remaining stock of dividents to households
    distribute_dividends_if!(
        indexfund_struct,
        all_hh,
        model
    )
end


function run_simulation(;
    T=400::Int,
    changed_params=nothing,
    full_output=true::Bool,
    labormarket_is_fordist=false::Bool,
    )::Tuple{Float64, Float64}

    to = TimerOutput()

    println("Init Model:")
    @timeit to "init" model, global_param, init_param, macro_struct, gov_struct, ep, labormarket_struct, indexfund_struct, climate_struct, cm_dat = initialize_model(T, labormarket_is_fordist; changed_params=changed_params)
    for t in 1:T

        if t % 100 == 0
            println("Step $t")
        end

        @timeit to "step" model_step!(
                                t,
                                to, 
                                global_param, 
                                init_param,
                                macro_struct, 
                                gov_struct,
                                ep,
                                labormarket_struct, 
                                indexfund_struct,
                                climate_struct,
                                cm_dat, 
                                model
                            )
    end

    @timeit to "save macro" save_macro_data(macro_struct)

    all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

    @timeit to "save findist" save_final_dist(all_hh, all_cp, all_kp, model)

    @timeit to "save climate" save_climate_data(ep, climate_struct, model)

    if full_output
        show(to)
        println()
    end

    return macro_struct.GDP[end], mean(macro_struct.GDP_growth)
end

run_simulation()