# import Pkg; Pkg.add("Distributions")
# import Pkg; Pkg.add("StatsBase") 
using Distributions
using StatsBase


mutable struct CapitalGoodProducer 
    id :: Int
    A_t :: Array{Float64}
    B_t :: Array{Float64}
    p_t :: Array{Float64}
    c_t :: Array{Float64}
    RD_t :: Array{Float64}
    IM_t :: Array{Float64}
    IN_t :: Array{Float64}
    S_t :: Array{Float64}
end


function innovate!(capital_good_producer, global_param, macro_struct, model_struct)
    """
    Checks if innovation is performed, then calls approapriate
        functions
    """

    set_RD(capital_good_producer, global_param, macro_struct)
    tech_choices = [(capital_good_producer.A_t[end], capital_good_producer.B_t[end])]

    # determine innovation of machines
    θ_IN = 1 - exp(-global_param.ζ * capital_good_producer.IN_t[end])
    if (Bernoulli(θ_IN) == 1)
        A_t_in = update_At(capital_good_producer, global_param)
        B_t_in = update_Bt(capital_good_producer, global_param)
        push!(tech_choices, (A_t_in, B_t_in))
    end

    # determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * capital_good_producer.IM_t[end])
    if (Bernoulli(θ_IM) == 1)
        A_t_im, B_t_im = imitate_technology(capital_good_producer, model_struct)
        push!(tech_choices, (A_t_im, B_t_im))
    end

    # make choice between possible technologies
    if length(tech_choices) == 1
        # if no new technologies, keep current technologies
        push!(capital_good_producer.A_t, capital_good_producer.A_t[end])
        push!(capital_good_producer.B_t, capital_good_producer.B_t[end])
    else
        c_h = map(x -> global_param.w[end]/x[0], tech_choices)
        p_h = map(x -> (1 + global_param.μ1)*x, c_h)
        r_h = p_h + global_param.b * c_h
        index = argmin(r_h)
        push!(capital_good_producer.A_t, tech_choices[index][0])
        push!(capital_good_producer.B_t, tech_choices[index][1])
    end

end


function imitate_technology(capital_good_producer, model_struct, F1)
    # Use inverse distances as weights for choice competitor,
    distances = model_struct.capital_good_euclidian_matrix[capital_good_producer.id, :]
    
    weights = zeros(F1)
    weights[1:capital_good_producer.id,:] = 1/distances[1:capital_good_producer.id]
    weights[capital_good_producer.id+1:end] = 1/distances[capital_good_producer.id+1:end]
    
    index = sample(1:F1, Weights(weights))

    A_t_im = model_struct.capital_good_producers[index].A_t[end]
    B_t_im = model_struct.capital_good_producers[index].B_t[end]
    
    return A_t_im, B_t_im
end


function update_At!(capital_good_producer, global_param)
    """
    Updates A_t
    """
    κ_A = Beta(global_param.α1, global_param.β1)
    A_t_in = A_t[end]*(1 + κ_A)
    return A_t_in
end


function update_Bt!(capital_good_producer, global_param) 
    κ_A = Beta(global_param.α1, global_param.β1)
    B_t_in = B_t[end]*(1 + κ_A)
    return B_t_in
end


function set_price(capital_good_producer, global_param, macro_struct)
    c_t = macro_struct.w[end] / capital_good_producer.B_t[end]
    p_t_new = (1 + global_param.μ1)c_t
    push!(capital_good_producer.p_t, p_t_new)
end


function set_RD(capital_good_producer, global_param, macro_struct)
    RD_t_new = global_param.ν*(capital_good_producer.S_t[end])
    push!(capital_good_producer.RD_t, RD_t_new)

    IN_t_new = global_param.ξ*RD_t_new
    IM_t_new = (1 - global_param.ξ)*RD_t_new
    push!(capital_good_producer.IN_t, IN_t_new)
    push!(capital_good_producer.IM_t, IN_t_new)
end