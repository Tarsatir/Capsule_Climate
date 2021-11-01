# import Pkg; Pkg.add("Distributions")
# import Pkg; Pkg.add("StatsBase") 
using Distributions
using StatsBase


mutable struct CapitalGoodProducer 
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
end


function innovate!(capital_good_producer, global_param, macro_struct, model_struct)
    """
    Checks if innovation is performed, then calls approapriate
        functions
    """

    set_RD(capital_good_producer, global_param, macro_struct)
    tech_choices = [(capital_good_producer.A[end], capital_good_producer.B[end])]

    # determine innovation of machines
    θ_IN = 1 - exp(-global_param.ζ * capital_good_producer.IN[end])
    if (Bernoulli(θ_IN) == 1)
        A_t_in = update_At(capital_good_producer, global_param)
        B_t_in = update_Bt(capital_good_producer, global_param)
        push!(tech_choices, (A_t_in, B_t_in))
    end

    # determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * capital_good_producer.IM[end])
    if (Bernoulli(θ_IM) == 1)
        A_t_im, B_t_im = imitate_technology(capital_good_producer, model_struct)
        push!(tech_choices, (A_t_im, B_t_im))
    end

    # make choice between possible technologies
    if length(tech_choices) == 1
        # if no new technologies, keep current technologies
        push!(capital_good_producer.A, capital_good_producer.A[end])
        push!(capital_good_producer.B, capital_good_producer.B[end])
    else
        c_h = map(x -> global_param.w[end]/x[0], tech_choices)
        p_h = map(x -> (1 + global_param.μ1)*x, c_h)
        r_h = p_h + global_param.b * c_h
        index = argmin(r_h)
        push!(capital_good_producer.A, tech_choices[index][0])
        push!(capital_good_producer.B, tech_choices[index][1])
    end

end


function send_brochures(capital_good_producer, model_struct, global_param)
    # TODO: reset brochures of consumer good producers each turn

    # set up brochure
    brochure = (capital_good_producer, capital_good_producer.p[end], capital_good_producer.c[end])

    # send brochure to historical clients
    for client in capital_good_producer.HC
        push!(client.brochures, brochure)
    end

    # select new clients, send brochure
    NC_potential = setdiff(model_struct.consumer_good_producers, consumer_good_producer.HC)
    indexes = rand(1:length(NC_potential), round(global_param.γ * length(consumer_good_producer.HC)))
    for index in indexes:
        push!(NC_potential[index].brochures, brochure)
    end 

end


function imitate_technology(capital_good_producer, model_struct, F1)
    # Use inverse distances as weights for choice competitor,
    distances = model_struct.capital_good_euclidian_matrix[capital_good_producer.id, :]
    
    weights = zeros(F1)
    weights[1:capital_good_producer.id,:] = 1/distances[1:capital_good_producer.id]
    weights[capital_good_producer.id+1:end] = 1/distances[capital_good_producer.id+1:end]
    
    index = sample(1:F1, Weights(weights))

    A_t_im = model_struct.capital_good_producers[index].A[end]
    B_t_im = model_struct.capital_good_producers[index].B[end]
    
    return A_t_im, B_t_im
end


function update_At!(capital_good_producer, global_param)
    """
    Updates A
    """
    κ_A = Beta(global_param.α1, global_param.β1)
    A_t_in = A[end]*(1 + κ_A)
    return A_t_in
end


function update_Bt!(capital_good_producer, global_param) 
    κ_A = Beta(global_param.α1, global_param.β1)
    B_t_in = B[end]*(1 + κ_A)
    return B_t_in
end


function set_price!(capital_good_producer, global_param, macro_struct)
    c_t = macro_struct.w[end] / capital_good_producer.B[end]
    p_t = (1 + global_param.μ1) * c_t
    push!(capital_good_producer.c, c_t)
    push!(capital_good_producer.p, p_t)
end


function set_RD!(capital_good_producer, global_param, macro_struct)
    RD_t_new = global_param.ν*(capital_good_producer.S[end])
    push!(capital_good_producer.RD, RD_t_new)

    IN_t_new = global_param.ξ*RD_t_new
    IM_t_new = (1 - global_param.ξ)*RD_t_new
    push!(capital_good_producer.IN, IN_t_new)
    push!(capital_good_producer.IM, IN_t_new)
end