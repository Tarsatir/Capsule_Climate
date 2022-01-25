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
include("labormarket.jl")
include("household.jl")
include("consumer_good.jl")
include("capital_good.jl")
include("macro.jl")



function initialize_model(;
    n_captlgood = 5,
    n_consrgood = 20,
    n_households = 250
    )

    # initialise model struct
    model = AgentBasedModel(Union{Household, CapitalGoodProducer, ConsumerGoodProducer})
    all_agents = All_Agents([], [], [], [])

    # initialize struct that holds global params
    global_param  = initialize_global_params()

    # initialize struct that holds macro variables
    macro_struct = initialize_macro()
    labormarket_struct = initialize_labormarket()

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

        push!(all_agents.households, hh)
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

        # initialize capital good stock
        machine_struct = initialize_machine()

        cp = initialize_cp(id, cp_id, machine_struct, n_consrgood)

        push!(all_agents.consumer_good_producers, cp)
        add_agent!(cp, model)
        id += 1
    end

    # initialize capital good producers
    for kp_id in 1:n_captlgood

        # make choice for historical clients
        HC = sample(all_agents.consumer_good_producers, 10; replace=false)

        # initialize capital good producer
        kp = initialize_kp(id, kp_id, HC, n_captlgood)

        push!(all_agents.capital_good_producers, kp)
        add_agent!(kp, model)
        id += 1
    end

    # determine distance matrix between capital good producers
    get_capgood_euclidian(all_agents, n_captlgood)

    # spread employed households over producers
    spread_employees_lm!(
        labormarket_struct, 
        all_agents.consumer_good_producers, 
        all_agents.capital_good_producers
    )

    return model, all_agents, global_param, macro_struct, labormarket_struct
end


function model_step!(model, all_agents, global_param, macro_struct, labormarket_struct)

    # reset brochures of all consumer good producers
    for cp in all_agents.consumer_good_producers
        reset_brochures_cp!(cp)
    end

    # (1) capital good producers innovate and send brochures
    for kp in all_agents.capital_good_producers
        innovate_kp!(kp, global_param, all_agents, macro_struct)
        send_brochures_kp!(kp, all_agents, global_param)
    end

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    for cp in all_agents.consumer_good_producers
        plan_production_cp!(cp, global_param)
        plan_investment_cp!(cp, global_param, all_agents.capital_good_producers)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp in all_agents.capital_good_producers
        plan_production_kp!(kp)
    end

    # (3) labor market matching process

    labormarket_process!(
        labormarket_struct, 
        all_agents.consumer_good_producers, 
        all_agents.capital_good_producers
    )

    

    # (4) producers hire workers. Government pays unemployment benefits
    # TODO

    # (5) Government receives income taxes. Households set consumption choice
    # TODO

    # (6) Households pick prefered products to buy and set budget and consumption package
    for hh in all_agents.households
        pick_cp_hh(hh, all_agents.consumer_good_producers)
    end

    # (6) Transactions take place on consumer market, consumer good producers
    # make up profits
    # TODO

    # (6) capital good producers deliver goods to consumer good producers
    # TODO

    # (7) government receives profit taxes
    # TODO

    # (7) macro-economic indicators are updated.
    # TODO
    # TODO update market share cp


end

to = TimerOutput()

@timeit to "init" model, all_agents, global_param, macro_struct, labormarket_struct = initialize_model()
for i in 1:1
    @timeit to "step" model_step!(model, all_agents, global_param, macro_struct, labormarket_struct)
end

show(to)
println()