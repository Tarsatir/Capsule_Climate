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
using BenchmarkTools
# using Setfield

include("model.jl")
include("misc.jl")
include("household.jl")
include("consumer_good.jl")
include("capital_good.jl")
include("macro.jl")



function initialize_model(;
    n_captlgood = 5,
    n_consrgood = 20,
    n_households = 100
    )

    # initialise model struct
    model = AgentBasedModel(Union{CapitalGoodProducer, ConsumerGoodProducer})
    all_agents = All_Agents([], [], [])

    # initialize struct that holds global params
    global_param  = initialize_global_params()

    # initialize struct that holds macro variables
    macro_struct = initialize_macro()

    # global id
    id = 0

    # initialize households
    for hh_id in 1:n_households
        
    end

    # initialize consumer good producers
    for cp_id in 1:n_consrgood

        # initialize capital good stock
        machine = Machine(
            1,                      # A: labor productivity machine
            0,                      # c: cost to produce machine
            40,                     # freq: freq machine owned by cp
            0                       # age: age of machine
        )

        consumer_good_producer = ConsumerGoodProducer(
            id,                     # global id
            cp_id,                  # cp id
            [],                     # p: hist prices
            [],                     # c: hist cost
            [],                     # RD: hist R&D spending
            [1000],                 # D: hist demand
            1000,                   # Dᵉ exp demand
            [100],                  # N: hist inventory
            200,                    # Nᵈ: desired inventory 
            [rand()],               # Q: hist production
            12000,                  # Qᵉ: exp production
            [rand()],               # I: hist investments
            [machine],              # Ξ: capital stock
            [25],                   # L: labor force
            25,                     # Lᵉ: exp labor force
            0,                      # ΔLᵈ: desired change in labor force
            [1.0],                  # w: wage level
            [],                     # brochures from kp
            [rand()],               # π: hist productivity
            [rand()],               # f: hist market share
            [0.05],                 # μ: hist markup
            [100000],               # Π: hist profits
            0,                      # cI: internal funds for investments
            0,                      # ΔDeb: changes in debt level
            Balance(               
                    0.0,            # - N: inventory
                    0.0,            # - K: capital
                    0.0,            # - NW: liquid assets
                    0.0,            # - Deb: debt
                    0.0             # - EQ: equity
                )
        )
        push!(all_agents.consumer_good_producers, consumer_good_producer)
        add_agent!(consumer_good_producer, model)
        id += 1
    end

    # initialize capital good producers
    for kp_id in 1:n_captlgood

        # make choice for historical clients
        HC = sample(all_agents.consumer_good_producers, 10; replace=false)

        # initialize capital good producer
        capital_good_producer = CapitalGoodProducer( # initial parameters based on rer98
            id,                     # global id
            kp_id,                  # kp id
            [1],                    # A: labor prod sold product
            [1],                    # B: labor prod own production
            [],                     # p: hist price data
            [],                     # c: hist cost data
            [],                     # RD: hist R&D expenditure
            [],                     # IM: hist immitation expenditure
            [],                     # IN: hist innovation expenditure
            [100],                  # S: hist revenue
            HC,                     # HC: hist clients
            [],                     # Π: hist profits
            [],                     # orders
            Balance(               
                    0.0,            # - N: inventory
                    0.0,            # - K: capital
                    1000.0,         # - NW: liquid assets
                    0.0,            # - Deb: debt
                    0.0             # - EQ: equity
                )               
        )
        push!(all_agents.capital_good_producers, capital_good_producer)
        add_agent!(capital_good_producer, model)
        id += 1
    end

    # determine distance matrix between capital good producers
    get_capgood_euclidian(all_agents, n_captlgood)

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
        plan_production_cp!(cp, global_param)
        plan_investment_cp!(cp, global_param)
    end

    # (2) capital good producers set labor demand based on ordered machines
    for kp in all_agents.capital_good_producers
        plan_production_kp!(kp, global_param)
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

    # (7) macro-economic indicators are updated.
    # TODO
    # TODO update market share cp


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

@time begin

model, all_agents, global_param, macro_struct = initialize_model()
model_step!(model, all_agents, global_param, macro_struct)

end