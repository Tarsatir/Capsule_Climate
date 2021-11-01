
mutable struct Machine
    A :: Float64
    B :: Float64
    c :: Float64
end

mutable struct ConsumerGoodProducer 
    id :: Int
    # A :: Array{Float64}
    # B :: Array{Float64}
    p :: Array{Float64}
    c :: Array{Float64}
    RD :: Array{Float64}
    D :: Array{Float64}
    N :: Array{Float64}
    K :: Array{Float64}
    I :: Array{Float64}
    Ksi :: Array
    brochures :: Array
end


function plan_investment!(consumer_good_producer, global_param, model_struct)

    # choose capital goods producer
    

    # compute investment
    EId = plan_expansion(consumer_good_producer, global_param, model_struct)
    RS = plan_replacement(consumer_good_producer, global_param, p_star, c_star)
    I_t = EId + RS
    push!(consumer_good_producer, EId, RS)
end


function plan_expansion(consumer_good_producer, global_param, model_struct)
    Nd = global_param.ι * consumer_good_producer.D[end]
    Qd = consumer_good_producer.D[end] + Nd - consumer_good_producer.N[end]
    Kd = Qd / consumer_good_producer.A[end]
    # TODO: check if this is the correct function
    EId = Kd - consumer_good_producer.K[end]

    return EId
end


function plan_replacement(consumer_good_producer, global_param, p_star, c_star)
    RS = Array{Float64}

    for machine in consumer_good_producer.Ksi
        if (p_star/(machine.A - c_star) <= global_param.b)
            push!(RS, machine)
        end
    end

    return RS
end