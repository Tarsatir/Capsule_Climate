# using Pkg; Pkg.add("Plots")

using Printf
using Statistics
using Distributions
using Random
using Agents

include("model.jl")
include("capital_good.jl")
include("consumer_good.jl")
include("macro.jl")


function initialize_model()
    model_struct = Model([], [], [])
    return model_struct
end

function initialize_global_params()
    global_param = GlobalParam(
        50,
        200,
        0.04,
        0.5,
        0.3,
        3.0,
        3.0,
        -0.15,
        0.15,
        0.5,
        0.04,
        0.1,
        3,
        20,
        0.04,
        1.0,
        1.0,
        2.0,
        0.01,
        0.1,
        0.9,
        0.1,
        0.9,
        2.0,
        4.0,
        1.0,
        0.0,
        0.0,
        0.1,
        0.4
    )
    return global_param
end

function initialize_macro()
    macro_struct = MacroEconomy([], [], [2], [], [])
    return macro_struct
end


model_struct = initialize_model()
global_param  = initialize_global_params()
macro_struct = initialize_macro()

# initialize capital good producers
for i in 1:global_param.F1
    capital_good_producer = CapitalGoodProducer(i, [rand()], [rand()], [], [], [], [], [], [rand()])
    push!(model_struct.capital_good_producers, capital_good_producer)
end

get_capgood_euclidian(model_struct, global_param)

for i in 1:global_param.F1
    capital_good_producer = model_struct.capital_good_producers[i]
    innovate!(capital_good_producer, global_param, macro_struct, model_struct)
end

# initialize consumer good producers
# for i in 1:global_param.F2
#     consumer_good_producer = ConsumerGoodProducer(i, [], [], [], [], [])
#     push!(model_struct.consumer_good_producers, consumer_good_producer)
# end

# println(model_struct.capital_good_producers)
# println(model_struct.consumer_good_producers)

