
mutable struct ConsumerGoodProducer 
    id :: Int
    A :: Array{Float64}
    B :: Array{Float64}
    p :: Array{Float64}
    c :: Array{Float64}
    RD :: Array{Float64}
    D :: Array{Float64}
    N :: Array{Float64}
    K :: Array{Float64}
end

function plan_production(consumer_good_producer, global_param, model_struct)
    Nd = global_param.ι * consumer_good_producer.D[end]
    Qd = consumer_good_producer.D[end] + Nd - consumer_good_producer.N[end]
    Kd = Qd / consumer_good_producer.A[end]
    # TODO: check if this is the correct function
    EId = Kd - consumer_good_producer.K[end]
end