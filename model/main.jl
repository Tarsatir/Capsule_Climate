# Install dependencies
# using Pkg; Pkg.add("Plots")
# import Pkg; Pkg.add("Distributions")
# import Pkg; Pkg.add("StatsBase") 
# import Pkg; Pkg.add("Agents")
# import Pkg; Pkg.add("Setfield")

using Printf
using Statistics
using Distributions
using StatsBase
using Random
using Agents
using BenchmarkTools
using TimerOutputs
# using Setfield

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
    n_captlgood = 5,
    n_consrgood = 20,
    n_households = 250
    )

    # initialise model struct
    model = AgentBasedModel(Union{Household, CapitalGoodProducer, ConsumerGoodProducer})
    all_agents = initialize_allagents()

    # initialize struct that holds global params
    global_param  = initialize_global_params()

    # initialize struct that holds macro variables
    macro_struct = initialize_macro()

    # initialize labor market struct
    labormarket_struct = initialize_labormarket()

    # initialize consumer market struct
    consumermarket_struct = initialize_consumermarket()

    # initialize government struct
    gov_struct = initialize_government()

    # global id
    id = 0

    # initialize households
    for hh_id in 1:n_households

        # determine if household will be employed
        employed = true
        if hh_id > 0.9 * n_households
            employed = false
        end

        hh = initialize_hh(id, hh_id, employed)

        push!(all_agents.all_hh, hh)
        add_agent!(hh, model)

        # add household to labor market struct based on employment status
        if hh.employed
            push!(labormarket_struct.employed, hh)
        else
            push!(labormarket_struct.unemployed, hh)
        end
    end

    # update unemployment rate
    update_unemploymentrate_lm(labormarket_struct)

    # initialize consumer good producers
    for cp_id in 1:n_consrgood

        # decide if producer makes basic or luxury goods
        type_good = "Basic"
        if cp_id > n_consrgood / 2
            type_good = "Luxury"
        end

        # initialize capital good stock
        machine_struct = initialize_machine()

        cp = initialize_cp(id, cp_id, machine_struct, n_consrgood, type_good)

        push!(all_agents.all_cp, cp)
        
        if type_good == "Basic"
            push!(all_agents.all_bp, cp)
        else
            push!(all_agents.all_lp, cp)
        end

        add_agent!(cp, model)
        id += 1
    end

    # initialize capital good producers
    for kp_id in 1:n_captlgood

        # make choice for historical clients
        HC = sample(all_agents.all_cp, 10; replace=false)

        # initialize capital good producer
        kp = initialize_kp(id, kp_id, HC, n_captlgood)

        push!(all_agents.all_kp, kp)
        add_agent!(kp, model)
        id += 1
    end

    # determine distance matrix between capital good producers
    get_capgood_euclidian(all_agents, n_captlgood)

    # spread employed households over producers
    spread_employees_lm!(
        labormarket_struct, 
        all_agents.all_cp, 
        all_agents.all_kp
    )

    return model, all_agents, global_param, macro_struct, gov_struct, labormarket_struct, consumermarket_struct
end


function model_step!( model, 
    all_agents, 
    global_param, 
    macro_struct, 
    gov_struct, 
    labormarket_struct,
    consumermarket_struct
    )

    # reset brochures of all consumer good producers
    for cp in all_agents.all_cp
        reset_brochures_cp!(cp)
    end

    # (1) capital good producers innovate and send brochures
    for kp in all_agents.all_kp
        innovate_kp!(kp, global_param, all_agents, macro_struct)
        send_brochures_kp!(kp, all_agents, global_param)
    end

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    for cp in all_agents.all_cp
        plan_production_cp!(cp, global_param)
        plan_investment_cp!(cp, global_param, all_agents.all_kp)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp in all_agents.all_kp
        plan_production_kp!(kp)
    end


    # (3) labor market matching process
    labormarket_process!(
        labormarket_struct, 
        all_agents.all_cp, 
        all_agents.all_kp,
        global_param.ϵ,
        gov_struct.UB
    )
    update_avg_T_unemp_lm(labormarket_struct)

    # (4) Producers pay workers their wage. Government pays unemployment benefits
    for p in vcat(all_agents.all_cp, all_agents.all_kp)
        pay_workers_p!(p)
    end

    pay_unemployment_benefits_gov!(gov_struct, labormarket_struct.unemployed)


    # (5) Production takes place for cp and kp
    for cp in all_agents.all_cp
        produce_goods_cp!(cp)
    end

    for kp in all_agents.all_kp
        produce_goods_kp!(kp)
    end


    # (5) Government receives income taxes
    levy_income_tax_gov!(gov_struct, all_agents.all_hh)
    # compute_budget_balance(gov_struct)

    # (6) Households pick prefered products to buy and set budget and consumption package
    for hh in all_agents.all_hh
        compute_exp_income_hh!(hh, 
                               labormarket_struct.P_HU, 
                               labormarket_struct.P_UU, 
                               gov_struct.UB)
        set_savingsrate_hh!(hh, labormarket_struct.avg_T_unemp, gov_struct.UB)
    end


    # (6) Transactions take place on consumer market
    consumermarket_process!(consumermarket_struct,
                            all_agents.all_hh,
                            all_agents.all_bp,
                            all_agents.all_lp,
                            gov_struct)

    # cp make up profits
    for cp in all_agents.all_cp
        compute_profits_cp!(cp)
    end

    # (6) kp deliver goods to cp, kp make up profits
    for kp in all_agents.all_kp
        send_orders_kp!(kp)
    end

    # (7) government receives profit taxes
    # TODO

    # (7) macro-economic indicators are updated.
    update_macro_stats(macro_struct, all_agents.all_hh, all_agents.all_cp, all_agents.all_kp)

    # TODO update market share cp


end

to = TimerOutput()

@timeit to "init" model, all_agents, global_param, macro_struct, gov_struct, labormarket_struct, consumermarket_struct = initialize_model()
for i in 1:50
    println("Step ", i)
    @timeit to "step" model_step!(model, all_agents, global_param, macro_struct, gov_struct, labormarket_struct, consumermarket_struct)
end

println(macro_struct.GDP)

show(to)
println()