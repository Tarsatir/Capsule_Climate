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
using Dates
# using UnPack


# Include files
include("../results/write_results.jl")
include("helpers/custom_schedulers.jl")
include("helpers/dist_matrix.jl")
include("helpers/tmpdata_storage.jl")
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
include("macro_markets/capitalmarket.jl")

include("helpers/modeldata_storage.jl")


"""
    initialize_model(T::Int64, labormarket_is_fordist::Bool, changed_params::Dict)

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
    T::Int64;
    changed_params::Union{Dict, Nothing},
    changedparams_ofat::Union{Dict, Nothing},
    changedtaxrates::Union{Vector, Nothing},
    track_firms_households::Bool=false
    )

    # Initialize scheduler
    scheduler = Agents.Schedulers.by_type(true, true)

    # Initialize struct that holds global params and initial parameters
    globalparam  = initialize_global_params(changed_params, changedparams_ofat)
    initparam = InitParam()

    # Initialize struct that holds macro variables
    macroeconomy = MacroEconomy(T=T)

    # Initialize labor market struct
    labormarket = LaborMarket()

    # Initialize government struct
    government = initgovernment(T, changedtaxrates)

    # Initialize energy producer
    ep = initialize_energy_producer(T, initparam, government.τᶜ, globalparam)

    # Initialize climate struct
    climate = Climate(T=T)


    # Initialize data structures for consumer- and capital market
    cmdata = CMData(
        n_hh=initparam.n_hh, 
        n_cp=initparam.n_cp
    )
    
    kmdata = KMData(
            zeros(Int64, initparam.n_kp, initparam.n_cp),
            zeros(Int64, initparam.n_kp, initparam.n_cp),
            # zeros(Int64, initparam.n_kp)  
        )

    properties = Dict(
                        :kp_brochures => Dict(),
                        :cmdata => cmdata,
                        :kmdata => kmdata,
                        :i_to_id => Dict("hh" => Dict(), "cp" => Dict(), "kp" => Dict())
                     )

    # Initialize model struct
    model = ABM(
                    Union{Household, CapitalGoodProducer, ConsumerGoodProducer};
                    properties=properties,
                    scheduler=scheduler, 
                    warn=false
                )

    # Global id
    id = 1

    # Initialize households
    hh_skills = sample_skills_hh(initparam)
    for hh_i in 1:initparam.n_hh

        hh = Household(id=id, skill=hh_skills[hh_i])
        hh.wʳ = max(government.w_min, hh.wʳ)
        add_agent!(hh, model)
        model.i_to_id["hh"][hh_i] = id

        id += 1
    end

    # Initialize consumer good producers
    for cp_i in 1:initparam.n_cp

        # Initialize capital good stock
        machines = initialize_machine_stock(
                        globalparam.freq_per_machine, 
                        initparam.n_machines_init,
                        η=globalparam.η; 
                        A_LP = initparam.A_LP_0,
                        A_EE = initparam.A_EE_0,
                        A_EF = initparam.A_EF_0
                   )

        t_next_update = 1
        if cp_i > 66
            t_next_update += 1
        end
        if cp_i > 132
            t_next_update += 1
        end

        cp = initialize_cp(
                id,
                cp_i,
                t_next_update,
                machines,  
                initparam.n_init_emp_cp, 
                globalparam.μ1,
                globalparam.ι,
                globalparam.b;
                n_consrgood=initparam.n_cp
            )
        update_n_machines_cp!(cp, globalparam.freq_per_machine)
        add_agent!(cp, model)
        model.i_to_id["cp"][cp_i] = id

        id += 1
    end

    # Initialize capital good producers
    for kp_i in 1:initparam.n_kp

        kp = initialize_kp(
                id, 
                kp_i, 
                initparam.n_kp,
                globalparam.b; 
                A_LP=initparam.A_LP_0,
                A_EE=initparam.A_EE_0,
                A_EF=initparam.A_EF_0, 
                B_LP=initparam.B_LP_0,
                B_EE=initparam.B_EE_0,
                B_EF=initparam.B_EF_0
             )
        add_agent!(kp, model)
        model.i_to_id["kp"][kp_i] = id

        # Initialize brochure of kp goods
        init_brochure!(kp, model)

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

    close_balance_all_p!(all_p, globalparam, government,
                         indexfund, 0, model)

    if track_firms_households
        firmdata = genfirmdata(all_cp, all_kp)
        householddata = genhouseholddata()
    else
        firmdata = nothing
        householddata = nothing
    end

    return model, globalparam, initparam, macroeconomy, government, ep, labormarket, 
           indexfund, climate, firmdata, householddata
end


function model_step!(
    t::Int64,
    t_warmup::Int64,
    to,
    globalparam::GlobalParam, 
    initparam::InitParam,
    macroeconomy::MacroEconomy, 
    government::Government,
    ep::EnergyProducer, 
    labormarket::LaborMarket,
    indexfund::IndexFund,
    climate::Climate,
    firmdata::Union{Nothing, DataFrame},
    householddata::Union{Nothing, Array},
    model::ABM
    )

    # Check if any global params are changed in ofat experiment
    check_changed_ofatparams(globalparam, t)

    # If end of warmup period reached, instate changed taxes
    if t == t_warmup
        instatetaxes!(government)
    end

    # Update schedulers
    @timeit to "schedule" all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

    # Redistrubute remaining stock of dividents to households
    @timeit to "distr div" distribute_dividends_if!(
        indexfund,
        government,
        all_hh,
        government.τᴷ,
        t,
        model
    )

    # Redistribute goverment balance
    resolve_gov_balance!(government, indexfund, globalparam, all_hh, model)

    # Update firm age
    for p_id in all_p
        model[p_id].age += 1
    end

    # Check if households still have enough bp and lp, otherwise sample more
    @timeit to "hh refill" for hh_id in all_hh
        refillsuppliers_hh!(model[hh_id], all_cp, initparam.n_cp_hh, model)
        resetincomes_hh!(model[hh_id])
    end

    # Clear current account, decide how many debts to repay, reset kp brochures of all cp
    @timeit to "clear account cp" for cp_id in all_cp
        model[cp_id].curracc = clear_firm_currentaccount_p!(model[cp_id].curracc)
        # reset_brochures_cp!(model[cp_id])
    end

    # (1) kp and ep innovate, and kp send brochures

    # Determine distance matrix between kp
    @timeit to "dist mat" kp_distance_matrix = get_capgood_euclidian(all_kp, model)

    @timeit to "kp innov" for kp_id in all_kp

        model[kp_id].curracc = clear_firm_currentaccount_p!(model[kp_id].curracc)

        innovate_kp!(
            model[kp_id],
            government, 
            globalparam, 
            all_kp, 
            kp_distance_matrix,
            macroeconomy.w̄_avg[max(t-1,1)],
            t,
            ep,
            model
        )

        # Update cost of production and price
        compute_c_kp!(model[kp_id], government, ep.pₑ[t])
        compute_p_kp!(model[kp_id])

        # Send brochures to cp
        send_brochures_kp!(model[kp_id], all_cp, globalparam, model)
    end

    # ep innovate, only when warmup period is left
    @timeit to "inn ep" innovate_ep!(ep, globalparam, t)

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    @timeit to "plan prod cp" for cp_id in all_cp
      
        # Plan production for this period
        # μ_avg = t > 1 ? macroeconomy.μ_cp[t-1] : globalparam.μ1
        plan_production_cp!(
            model[cp_id],
            government, 
            ep,
            globalparam, 
            government.τˢ,
            length(all_hh),
            length(all_cp),
            t, 
            model
        )

        if model[cp_id].t_next_update == t  
            model[cp_id].t_next_update += globalparam.update_period
        end


        # Reset desired and ordered machines
        reset_desired_ordered_machines_cp!(model[cp_id])

        # Plan investments for this period
        # plan_investment_cp!(model[cp_id], government, all_kp, globalparam, ep, t, model)

        # Rank producers based on cost of acquiring machines
        rank_producers_cp!(model[cp_id], government, globalparam.b, all_kp, ep, t, model)

        # Update expected long-term production
        update_Qᵉ_cp!(model[cp_id], globalparam.ω, globalparam.ι)
    end

    # (2) capital good producers set labor demand based on expected ordered machines
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

    # Update kp production capacitity
    for kp_id in all_kp
        update_prod_cap_kp!(model[kp_id], globalparam)
    end

    for cp_id in all_cp
        check_funding_restrictions_cp!(model[cp_id], government, globalparam, ep.pₑ[t])
    end

    # Let cp order capital goods from kp
    @timeit to "capitalmarket" capitalmarket_process!(
                all_cp, 
                all_kp,
                government,
                globalparam,
                ep,
                t, 
                model
            )


    # (4) Producers pay workers their wage. Government pays unemployment benefits
    @timeit to "pay workers" for p_id in all_p
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
    produce_energy_ep!(
        ep, 
        government, 
        all_cp, 
        all_kp, 
        globalparam, 
        indexfund,
        initparam.frac_green, 
        t, 
        t_warmup,
        model
    )
    
    # (6) Transactions take place on consumer market

    all_W = map(hh_id -> model[hh_id].W, all_hh)

    # TODO: change to actual expected returns
    ERt = t == 1 ? 0.07 : macroeconomy.returns_investments[t-1]

    # Households set consumption budget
    Wmin = minimum(all_W)
    Wmax = maximum(all_W)
    Wmedian = median(all_W)
    @timeit to "set budget" @inbounds for hh_id in all_hh
        set_consumption_budget_hh!(
            model[hh_id], 
            government.UB, 
            globalparam,
            ERt,
            labormarket.P_getunemployed,
            labormarket.P_getemployed,
            Wmin,
            Wmax,
            Wmedian,
            model
        )
    end

    # Consumer market process
    @timeit to "consumermarket" consumermarket_process!(
        all_hh,
        all_cp,
        government,
        globalparam,
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

    # println("   $t 7 start $(Dates.format(now(), "HH:MM"))")

    # (6) kp deliver goods to cp, kp make up profits
    @timeit to "send machines kp" for kp_id in all_kp
        send_ordered_machines_kp!(model[kp_id], ep, globalparam, t, model)
    end

    # Close balances of all firms
    for cp_id in all_cp
        # Update amount of owned capital, increase machine age
        update_n_machines_cp!(model[cp_id], globalparam.freq_per_machine)
        increase_machine_age_cp!(model[cp_id])
    end 

    # Close balances of firms, if insolvent, liquidate firms
    @timeit to "close balance" close_balance_all_p!(
                                all_p, 
                                globalparam,
                                government,
                                indexfund, 
                                t, 
                                model
                            )

    # println("   $t 7 start $(Dates.format(now(), "HH:MM"))")

    # (7) government receives profit taxes and computes budget balance
    levy_profit_tax_gov!(government, all_p, t, model)
    compute_budget_balance(government, t)

    # Update market shares of cp and kp
    update_marketshare_p!(all_cp, model)
    update_marketshare_p!(all_kp, model)
    
    # Select producers that will be declared bankrupt and removed
    @timeit to "check br" bankrupt_cp, bankrupt_kp, bankrupt_kp_i = check_bankrupty_all_p!(all_p, all_kp, globalparam, model)

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

    if firmdata ≠ nothing
        appendfirmdata!(firmdata, all_cp, all_kp, t, model)
    end

    if householddata ≠ nothing
        appendhouseholddata!(householddata, all_hh, t, t_warmup, model)
    end

    # Update climate parameters, compute new carbon equilibria and temperature change
    collect_emissions_cl!(climate, all_cp, all_kp, ep, t, globalparam.t_warmup, model)

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

    return model, globalparam, initparam, macroeconomy, government, ep, labormarket, 
            indexfund, climate, ep, firmdata, householddata
end


"""
    run_simulation(T::Int64, changed_params::Bool, full_output::Bool)

## Performs a full simulation.
    - Initializes model and agent structs.
    - Runs model `T` time steps.
    - Writes simulation results to csv.
"""
function run_simulation(;
    T::Int64=660,
    t_warmup::Int64=300,
    changed_params::Union{Dict,Nothing}=nothing,
    changedparams_ofat::Union{Dict,Nothing}=nothing,
    changedtaxrates::Union{Vector,Nothing}=nothing,
    full_output::Bool=true,
    threadnr::Int64=1,
    sim_nr::Int64=0,
    savedata::Bool=false,
    track_firms_households::Bool=false,
    seed::Int64=Random.rand(1000:9999)
    )

    # Set seed of simulation
    Random.seed!(seed)

    println("thread $(Threads.threadid()), sim $sim_nr has started on $(Dates.format(now(), "HH:MM"))")

    to = TimerOutput()

    @timeit to "init" model, globalparam, initparam, macroeconomy, government, ep, labormarket, 
                      indexfund, climate, firmdata, householddata = initialize_model(
                            T; 
                            changed_params=changed_params,
                            changedparams_ofat=changedparams_ofat, 
                            changedtaxrates=changedtaxrates, 
                            track_firms_households=track_firms_households
                    )
    
    @time for t in 1:T
        @timeit to "step" model_step!(
                                t,
                                t_warmup,
                                to, 
                                globalparam, 
                                initparam,
                                macroeconomy, 
                                government,
                                ep,
                                labormarket, 
                                indexfund,
                                climate,
                                # cmdata,
                                firmdata,
                                householddata, 
                                model
                            )
    end

    if savedata
        @timeit to "save macro" save_macro_data(macroeconomy)

        all_hh, all_cp, all_kp, all_p = schedule_per_type(model)

        @timeit to "save findist" save_final_dist(all_hh, all_cp, all_kp, model)

        @timeit to "save climate" save_climate_data(ep, climate, model)

        if track_firms_households
            save_household_quartiles(householddata)
        end
    end

    if full_output
        show(to)
        println()
    end

    # Pack output in struct

    println("thread $threadnr, sim $sim_nr has finished on $(Dates.format(now(), "HH:MM"))")

    if !track_firms_households
        return genrunoutput(macroeconomy, ep, climate)
    else
        return genrunoutput(macroeconomy, ep, climate), firmdata, householddata
    end

end

@time run_simulation(
    savedata=true,
    track_firms_households=true,
    seed=1234    
)

nothing