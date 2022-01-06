# Install dependencies
# using Pkg; Pkg.add("Plots")
# import Pkg; Pkg.add("Distributions")
# import Pkg; Pkg.add("StatsBase") 
# import Pkg; Pkg.add("Agents")
# import Pkg; Pkg.add("Setfield")

using Printf
using Statistics
using Distributions
using Random
using Agents
# using Setfield

include("model.jl")
include("capital_good.jl")
include("consumer_good.jl")
include("macro.jl")


function initialize_model(;
    n_captlgood = 5,
    n_consrgood = 20,
    n_households = 1000
    )

    # initialise model struct
    model = AgentBasedModel(Union{CapitalGoodProducer, ConsumerGoodProducer})
    all_agents = All_Agents([], [], [])

    # initialize struct that holds global params
    global_param  = initialize_global_params()

    # initialize struct that holds macro variables
    macro_struct = initialize_macro()

    # initialize capital good producers
    id = 0
    for kp_id in 1:n_captlgood
        capital_good_producer = CapitalGoodProducer(
            id,
            kp_id, 
            [rand()], 
            [rand()], 
            [rand()], 
            [rand()], 
            [], 
            [], 
            [], 
            [100], 
            [],
            [],
            []
        )
        push!(all_agents.capital_good_producers, capital_good_producer)
        add_agent!(capital_good_producer, model)
        id += 1
    end

    # determine distance matrix between capital good producers
    get_capgood_euclidian(all_agents, n_captlgood)

    # initialize consumer good producers
    for _ in 1:n_consrgood
        consumer_good_producer = ConsumerGoodProducer(
            id, 
            [], 
            [], 
            [], 
            [], 
            [rand()], 
            [rand()], 
            [], 
            [], 
            [rand()],
            [rand()],
            [100000],
            0,
            0,
            Balance(0, 0, 0, 0, 0),
        )
        push!(all_agents.consumer_good_producers, consumer_good_producer)
        add_agent!(consumer_good_producer, model)
        id += 1
    end

    return model, all_agents, global_param, macro_struct
end


function model_step!(model, all_agents, global_param, macro_struct)

    # reset brochures of all consumer good producers
    for cp in all_agents.consumer_good_producers
        reset_brochures!(cp)
    end

    # (1) capital good producers innovate and send brochures
    for kp in all_agents.capital_good_producers
        innovate!(kp, global_param, all_agents, macro_struct)
        send_brochures!(kp, all_agents, global_param)
    end

    # (2) consumer good producers estimate demand, set production and set
    # demand for L and K
    for cp in all_agents.consumer_good_producers
        println(cp.brochures)
        # TODO
    end

    # (3) workers apply for jobs
    # TODO

    # (4) producers hire workers. Government pays unemployment benefits
    # TODO

    # (5) Government receives income taxes. Households set consumption choice
    # TODO

    # (6) Transactions take place on consumer market, consumer good producers
    # make up profits
    # TODO

    # (6) capital good producers deliver goods to consumer good producers
    # TODO

    # (7) government receives profit taxes
    # TODO

    # (7) macro-economic indicators are updated
    # TODO


end


# function model_step!(model_struct)

    # # reset brochures of all consumer good producers
    # for consumer_good_producer in model_struct.consumer_good_producers
    #     consumer_good_producer.brochures = []
    # end

    # # let capital good producers innovate and send brochures
    # for capital_good_producer in model_struct.capital_good_producers
    #     innovate!(capital_good_producer, global_param, macro_struct, model_struct)
    #     send_brochures!(capital_good_producer, model_struct, global_param)
    # end

#     for consumer_good_producer in model_struct.consumer_good_producers
#         # let consumer good producers set investments and send orders
#         plan_investment!(consumer_good_producer, global_param, model_struct)

#         # let consumer good producers determine prices and competetiveness
#         w_t = macro_struct.w[end]
#         l_t = macro_struct.l[end]
#         compute_p_c_E!(consumer_good_producer, global_params, w_t, l_t)
#     end

#     # compute the average competetiveness of all consumer good producers
#     E_bar = compute_average_competitiveness!(macro_struct, model.struct.consumer_good_producers)

#     # let consumer good producers update market share, true demand and profits
#     for consumer_good_producer in model_struct.consumer_good_producers
#         set_production!(consumer_good_producer, E_bar, global_param.χ, 
#                         macro_struct.C[end], macro_struct.r)
#     end
# end

model, all_agents, global_param, macro_struct = initialize_model()

model_step!(model, all_agents, global_param, macro_struct)
