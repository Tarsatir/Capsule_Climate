using Statistics
using Distributions
using StatsBase
using Random
using Agents
using BenchmarkTools
using TimerOutputs
using DataStructures
using Parameters
using SparseArrays
using PyCall

# using Conda
# Conda.add("scipy")

# scipy = pyimport("scipy")

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

include("macro_markets/labormarket.jl")
include("macro_markets/consumermarket.jl")


"""
    initialize_model(T::Int, labormarket_is_fordist::Bool, changed_params::Dict)

## Initializes model.
Receives:
* `T`, total number of time steps.
* `changed_params`, changed global parameters in case
    of run with other parameters than default parameters.
    
Returns:
* `model`: Agent Based Model struct
* `globalparam`: struct containing global parameter values
* `macroeconomy`: mutable struct containing macro variables
* `government`: mutable struct containing variables concerning the government
* `labormarket`: mutable struct containing variables concerning the labor market
"""
function initialize_model(
    T::Int;
    changed_params::Union{Dict, Nothing},
    changedtaxrates::Union{Vector, Nothing}
    )

    # Initialise model struct
    scheduler = Agents.Schedulers.by_type(true, true)
    model = ABM(Union{Household, CapitalGoodProducer, ConsumerGoodProducer};
                scheduler=scheduler, warn=false)

    # Initialize struct that holds global params and initial parameters
    globalparam  = initialize_global_params(changed_params)
    initparam = InitParam()

    # Initialize struct that holds macro variables
    macroeconomy = MacroEconomy(T=T)

    # Initialize labor market struct
    labormarket = LaborMarket()

    # Initialize government struct
    government = initgovernment(T, changedtaxrates)

    # Initialize energy producer
    ep = initialize_energy_producer(T, initparam, globalparam)

    # Initialize climate struct
    climate = Climate(T=T)

    # Global id
    id = 1

    # Initialize households
    hh_skills = sample_skills_hh(initparam)
    for i in 1:initparam.n_hh

        hh = Household(id=id, skill=hh_skills[i])
        hh.wʳ = max(government.w_min, hh.wʳ)
        add_agent!(hh, model)

        id += 1
    end

    # Initialize consumer good producers
    for i in 1:initparam.n_cp

        # Initialize capital good stock
        machines = initialize_machine_stock(globalparam.freq_per_machine, 
                        initparam.n_machines_init,
                        η=globalparam.η; 
                        A_LP = initparam.A_LP_0,
                        A_EE = initparam.A_EE_0,
                        A_EF = initparam.A_EF_0
                   )

        t_next_update = 1
        if i > 66
            t_next_update += 1
        end
        if i > 132
            t_next_update += 1
        end

        cp = initialize_cp(
                id,
                t_next_update,
                machines,  
                initparam.n_init_emp_cp, 
                globalparam.μ1,
                globalparam.ι;
                n_consrgood=initparam.n_cp
            )
        update_n_machines_cp!(cp, globalparam.freq_per_machine)
        add_agent!(cp, model)

        id += 1
    end

    # Initialize capital good producers
    for kp_i in 1:initparam.n_kp

        kp = initialize_kp(
                id, 
                kp_i, 
                initparam.n_kp; 
                A_LP=initparam.A_LP_0,
                A_EE=initparam.A_EE_0,
                A_EF=initparam.A_EF_0, 
                B_LP=initparam.B_LP_0,
                B_EE=initparam.B_EE_0,
                B_EF=initparam.B_EF_0
             )
        add_agent!(kp, model)

        id += 1
    end

    # Initialize schedulers
    all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

    # Initialize index fund struct
    indexfund = IndexFund(T=T)

    # Let all capital good producers select historical clients
    for kp_id in all_kp
        select_HC_kp!(model[kp_id], all_cp)
    end

    # Let all households select cp of both types for trading network
    for hh_id in all_hh
        select_cp_hh!(model[hh_id], all_cp, initparam.n_cp_hh)
    end

    # Spread employed households over producerss
    spread_employees_lm!(
        labormarket,
        government,
        all_hh, 
        all_cp, 
        all_kp,
        initparam.n_init_emp_cp,
        initparam.n_init_emp_kp,
        model
    )

    for p_id in all_p
        update_mean_skill_p!(model[p_id], model)
    end

    close_balance_all_p!(all_p, globalparam, government.τᴾ,
                         indexfund, 0, model)

    cm_dat = CMData(n_hh=initparam.n_hh, n_cp=initparam.n_cp)

    return model, globalparam, initparam, macroeconomy, government, ep, labormarket, indexfund, climate, cm_dat
end


function model_step!(
    t::Int,
    to,
    globalparam::GlobalParam, 
    initparam::InitParam,
    macroeconomy::MacroEconomy, 
    government::Government,
    ep::EnergyProducer, 
    labormarket::LaborMarket,
    indexfund::IndexFund,
    climate::Climate,
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
    for hh_id in all_hh
        refillsuppliers_hh!(model[hh_id], all_cp, initparam.n_cp_hh, model)
        resetincomes_hh!(model[hh_id])
    end

    # Clear current account, decide how many debts to repay, reset kp brochures of all cp
    for cp_id in all_cp
        model[cp_id].curracc = clear_firm_currentaccount_p!(model[cp_id].curracc)
        reset_brochures_cp!(model[cp_id])
    end

    # (1) kp and ep innovate, and kp send brochures

    # Determine distance matrix between kp
    @timeit to "dist mat" kp_distance_matrix = get_capgood_euclidian(all_kp, model)

    @timeit to "kp innov" for kp_id in all_kp

        model[kp_id].curracc = clear_firm_currentaccount_p!(model[kp_id].curracc)

        innovate_kp!(
            model[kp_id], 
            globalparam, 
            all_kp, 
            kp_distance_matrix,
            macroeconomy.w̄_avg[max(t-1,1)],
            t,
            ep,
            model
        )

        # Update cost of production and price
        compute_c_kp!(model[kp_id])
        compute_p_kp!(model[kp_id])

        # Send brochures to cp
        send_brochures_kp!(model[kp_id], all_cp, globalparam, model)
    end

    # ep innovate
    innovate_ep!(ep, globalparam, t)

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    @timeit to "plan prod inv cp" for cp_id in all_cp
      
        # Plan production for this period
        μ_avg = t > 1 ? macroeconomy.μ_cp[t-1] : globalparam.μ1
        plan_production_cp!(
            model[cp_id], 
            globalparam, 
            μ_avg, 
            government.τˢ, 
            t, 
            model
        )

        if model[cp_id].t_next_update == t  
            model[cp_id].t_next_update += globalparam.update_period
        end

        # Plan investments for this period
        plan_investment_cp!(model[cp_id], all_kp, globalparam, ep, t, model)
    end

    # (2) capital good producers set labor demand based on ordered machines
    @timeit to "plan prod kp"  for kp_id in all_kp
        plan_production_kp!(model[kp_id], globalparam, model)
    end

    # (3) labor market matching process
    @timeit to "labormarket" labormarket_process!(
        labormarket,
        all_hh, 
        all_p,
        globalparam,
        government,
        t, 
        model,
        to
    )

    # Update mean skill level of employees
    for p_id in all_p
        update_mean_skill_p!(model[p_id], model)
    end


    # (4) Producers pay workers their wage. Government pays unemployment benefits
    for p_id in all_p
        pay_workers_p!(
            model[p_id],
            government,
            t, 
            model
        )
    end
    pay_unemployment_benefits_gov!(government, labormarket.unemployed, t, model)


    # (5) Production takes place for cp and kp

    # Update energy prices
    compute_pₑ_ep!(ep, t)

    @timeit to "consumer prod" for cp_id in all_cp
        produce_goods_cp!(model[cp_id], ep, globalparam, t)
    end

    @timeit to "capital prod" for kp_id in all_kp
        produce_goods_kp!(model[kp_id], ep, globalparam, t)
    end

    # Let energy producer meet energy demand
    produce_energy_ep!(ep, all_cp, all_kp, globalparam, indexfund, t, model)

    
    # (6) Transactions take place on consumer market

    all_I = map(hh_id -> model[hh_id].total_I, all_hh)

    # Households set consumption budget
    @timeit to "set budget" @inbounds for hh_id in all_hh
        set_consumption_budget_hh!(model[hh_id], all_I, globalparam, model)
    end

    # Consumer market process
    @timeit to "consumermarket" consumermarket_process!(
        all_hh,
        all_cp,
        government,
        globalparam,
        cm_dat,
        t,
        model,
        to
    )

    # Households decide to switch suppliers based on satisfied demand and prices
    @timeit to "decide switching hh" decide_switching_all_hh!(
        globalparam,
        all_hh,
        all_cp,
        all_p,
        initparam.n_cp_hh,
        model,
        to
    )

    # (6) kp deliver goods to cp, kp make up profits
    for kp_id in all_kp
        send_ordered_machines_kp!(model[kp_id], ep, globalparam, t, model)
    end

    # Close balances of all firms
    for cp_id in all_cp
        # Update amount of owned capital, increase machine age
        update_n_machines_cp!(model[cp_id], globalparam.freq_per_machine)
        increase_machine_age_cp!(model[cp_id])
    end 

    # Close balances of firms, if insolvent, liquidate firms
    @timeit to "close balance" close_balance_all_p!(all_p, globalparam, government.τᴾ,
                                    indexfund, t, model)

    # (7) government receives profit taxes and computes budget balance
    levy_profit_tax_gov!(government, all_p, t, model)
    compute_budget_balance(government, t)
    resolve_gov_balance!(government, indexfund, all_hh, model)

    # Update market shares of cp and kp
    update_marketshare_p!(all_cp, model)
    update_marketshare_p!(all_kp, model)
    
    # Select producers that will be declared bankrupt and removed
    bankrupt_cp, bankrupt_kp, bankrupt_kp_i = check_bankrupty_all_p!(all_p, all_kp, globalparam, model)

    # (7) macro-economic indicators are updated.
    @timeit to "update macro ts" update_macro_timeseries(
        macroeconomy,
        t, 
        all_hh, 
        all_cp, 
        all_kp,
        all_p,
        ep,
        bankrupt_cp,
        bankrupt_kp,
        labormarket,
        government,
        indexfund,
        globalparam,
        model,
        to
    )

    # Update climate parameters, compute new carbon equilibria and temperature change
    collect_emissions_cl!(climate, all_cp, all_kp, ep, t, model)

    # NOTE: climate process no longer tracked
    # carbon_equilibrium_tempchange_cl!(climate, t)

    # Remove bankrupt companies.
    @timeit to "kill bankr p" kill_all_bankrupt_p!(
        bankrupt_cp,
        bankrupt_kp, 
        all_hh,
        all_kp,
        labormarket,
        indexfund, 
        model
    )
    update_unemploymentrate_lm!(labormarket)


    # Replace bankrupt kp. Do this before you replace cp, such that new cp can also choose
    # from new kp
    @timeit to "replace kp" replace_bankrupt_kp!(
        bankrupt_kp, 
        bankrupt_kp_i, 
        all_kp,
        globalparam,
        indexfund,
        initparam,
        macroeconomy,
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
        globalparam,
        indexfund,
        macroeconomy,
        t,
        model
    )

    # Redistrubute remaining stock of dividents to households
    distribute_dividends_if!(
        indexfund,
        government,
        all_hh,
        government.τᴷ,
        t,
        model
    )
end


"""
    run_simulation(T::Int, changed_params::Bool, full_output::Bool)

## Performs a full simulation.
    - Initializes model and agent structs.
    - Runs model `T` time steps.
    - Writes simulation results to csv.
"""
function run_simulation(;
    T=460::Int,
    changed_params=nothing,
    changedtaxrates::Union{Vector,Nothing}=nothing,
    full_output::Bool=true,
    threadnr::Int64=1,
    savedata::Bool=false
    )

    to = TimerOutput()

    @timeit to "init" model, globalparam, initparam, macroeconomy, government, ep, labormarket, indexfund, climate, cm_dat = initialize_model(T; changed_params=changed_params, changedtaxrates=changedtaxrates)
    for t in 1:T

        @timeit to "step" model_step!(
                                t,
                                to, 
                                globalparam, 
                                initparam,
                                macroeconomy, 
                                government,
                                ep,
                                labormarket, 
                                indexfund,
                                climate,
                                cm_dat, 
                                model
                            )
    end

    if savedata
        @timeit to "save macro" save_macro_data(macroeconomy)

        all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

        @timeit to "save findist" save_final_dist(all_hh, all_cp, all_kp, model)

        @timeit to "save climate" save_climate_data(ep, climate, model)
    end

    if full_output
        show(to)
        println()
    end

    # Pack output in struct
    @timeit to "save output" runoutput = RunOutput(
        macroeconomy.GDP_growth,
        macroeconomy.U,
        macroeconomy.GINI_I,
        macroeconomy.GINI_W,
        macroeconomy.FGT,
        macroeconomy.avg_π_LP,
        macroeconomy.avg_π_EE,
        macroeconomy.avg_π_EF,
        climate.emissions_index
    )
    return runoutput
end

# run_simulation(savedata=true)