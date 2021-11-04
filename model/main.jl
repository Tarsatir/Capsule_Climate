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


function initialize()

    # initialise model, global_param and macro_structs
    model_struct = initialize_model()
    global_param  = initialize_global_params()
    macro_struct = initialize_macro()

    # initialize capital good producers
    for i in 1:global_param.F1
        capital_good_producer = CapitalGoodProducer(i, [rand()], [rand()], [rand()], [rand()], [], [], [], [rand()], [])
        push!(model_struct.capital_good_producers, capital_good_producer)
    end

    # determine distance matrix between capital good producers
    get_capgood_euclidian(model_struct, global_param)

    # initialize consumer good producers
    for i in 1:global_param.F2
        consumer_good_producer = ConsumerGoodProducer(i, [], [], [], [], [rand()], [rand()], [], [], [])
        push!(model_struct.consumer_good_producers, consumer_good_producer)
    end

    return model_struct, global_param, macro_struct
end


function iterate(model_struct)

    # reset brochures of all consumer good producers
    for consumer_good_producer in model_struct.consumer_good_producers
        consumer_good_producer.brochures = []
    end

    # let capital good producers innovate and send brochures
    for capital_good_producer in model_struct.capital_good_producers
        innovate!(capital_good_producer, global_param, macro_struct, model_struct)
        send_brochures!(capital_good_producer, model_struct, global_param)
    end

    for consumer_good_producer in model_struct.consumer_good_producers
        # let consumer good producers set investments and send orders
        plan_investment!(consumer_good_producer, global_param, model_struct)

        # let consumer good producers determine prices and competetiveness
        w_t = macro_struct.w[end]
        l_t = macro_struct.l[end]
        compute_p_c_E!(consumer_good_producer, global_params, w_t, l_t)
    end

    # compute the average competetiveness of all consumer good producers
    E_bar = compute_average_competitiveness!(macro_struct, model.struct.consumer_good_producers)

    # let consumer good producers update market share, true demand and profits
    for consumer_good_producer in model_struct.consumer_good_producers
        set_production!(consumer_good_producer, E_bar, global_param.χ, 
                        macro_struct.C[end], macro_struct.r)
    end


end

model_struct, global_param, macro_struct = initialize()

# println(model_struct)

# println(model_struct.capital_good_producers)
# println(model_struct.consumer_good_producers)

iterate(model_struct)
