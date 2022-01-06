using Distributions
using StatsBase


mutable struct Machine
    A :: Float64
    c :: Float64
    freq :: Int16
    age :: Float16
end


mutable struct CapitalGoodProducer <: AbstractAgent
    id :: Int
    A :: Array{Float64}
    B :: Array{Float64}
    p :: Array{Float64}
    c :: Array{Float64}
    RD :: Array{Float64}
    IM :: Array{Float64}
    IN :: Array{Float64}
    S :: Array{Float64}
    HC :: Array{Float64}
    Π :: Array{Float64}
    orders :: Array
end


function innovate!(;
    kp, 
    global_param, 
    macro_struct, 
    model_struct
    )
    """
    Checks if innovation is performed, then calls approapriate
        functions
    """
    
    # determine levels of R&D, and how to divide under IN and IM
    set_RD!(kp, global_param)
    tech_choices = [(kp.A[end], kp.B[end])]

    # determine innovation of machines (Dosi et al (2010); eq. 4)
    θ_IN = 1 - exp(-global_param.ζ * kp.IN[end])
    if (Bernoulli(θ_IN) == 1)
        A_t_in = update_At(global_param)
        B_t_in = update_Bt(global_param)
        push!(tech_choices, (A_t_in, B_t_in))
    end

    # determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * kp.IM[end])
    if (Bernoulli(θ_IM) == 1)
        A_t_im, B_t_im = imitate_technology(kp, model_struct)
        push!(tech_choices, (A_t_im, B_t_im))
    end

    # make choice between possible technologies
    if length(tech_choices) == 1
        # if no new technologies, keep current technologies
        push!(kp.A, kp.A[end])
        push!(kp.B, kp.B[end])
    else
        c_h = map(x -> global_param.w[end]/x[0], tech_choices)
        p_h = map(x -> (1 + global_param.μ1)*x, c_h)
        r_h = p_h + global_param.b * c_h
        index = argmin(r_h)
        push!(kp.A, tech_choices[index][0])
        push!(kp.B, tech_choices[index][1])
        push!(kp.p, p_h[index])
        push!(kp.c, c_h[index])
    end

end


function send_brochures!(kp, model_struct, global_param)

    # set up brochure
    brochure = (kp, kp.p[end], kp.c[end])

    # send brochure to historical clients
    for client in kp.HC
        push!(client.brochures, brochure)
    end

    # select new clients, send brochure
    NC_potential = setdiff(model_struct.consumer_good_producers, kp.HC)

    n_choices = round(global_param.γ * length(kp.HC))
    if length(kp.HC) == 0
        n_choices = 10
    end
    
    NC = sample(NC_potential, n_choices, replace=false)
    for consumer_good_producer in NC
        push!(consumer_good_producer.brochures, brochure)
    end 

end


function imitate_technology(kp, model_struct, F1)
    # use inverse distances as weights for choice competitor,
    distances = model_struct.capital_good_euclidian_matrix[kp.id, :]
    
    weights = zeros(F1)
    weights[1:kp.id,:] = 1/distances[1:kp.id]
    weights[kp.id+1:end] = 1/distances[kp.id+1:end]
    
    index = sample(1:F1, Weights(weights))

    A_t_im = model_struct.capital_good_producers[index].A[end]
    B_t_im = model_struct.capital_good_producers[index].B[end]
    
    return A_t_im, B_t_im
end


function set_production!(kp)
    production = 0
    for order in kp.orders
        production += order[2]
    end
end


function update_At!(global_param)
    # determines new labor productivity of machine produced for cp
    κ_A = Beta(global_param.α1, global_param.β1)
    A_t_in = A[end]*(1 + κ_A)
    return A_t_in
end


function update_Bt!(global_param)
    # determines new labor productivity of own production method 
    κ_A = Beta(global_param.α1, global_param.β1)
    B_t_in = B[end]*(1 + κ_A)
    return B_t_in
end


function set_price!(kp, global_param, macro_struct)
    c_t = macro_struct.w[end] / kp.B[end]
    p_t = (1 + global_param.μ1) * c_t
    push!(kp.c, c_t)
    push!(kp.p, p_t)
end


"""
Determines the level of R&D, and how it is divided under innovation (IN) 
and immitation (IM). based on Dosi et al (2010)
"""
function set_RD!(kp, global_param)
    # determine total R&D spending at time t, (Dosi et al, 2010; eq. 3)
    RD_new = global_param.ν*(kp.S[end])
    push!(kp.RD, RD_new)

    # decide fractions innovation (IN) and immitation (IM), 
    # (Dosi et al, 2010; eq. 3.5)
    IN_t_new = global_param.ξ*RD_new
    IM_t_new = (1 - global_param.ξ)*RD_new
    push!(kp.IN, IN_t_new)
    push!(kp.IM, IM_t_new)
end


function send_order(consumer_good_producer, order)

end