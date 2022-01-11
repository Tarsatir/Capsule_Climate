using Distributions
using StatsBase


mutable struct CapitalGoodProducer <: AbstractAgent
    id :: Int                   # global id
    kp_id :: Int                # kp id
    A :: Array{Float64}         # labor prod sold product
    B :: Array{Float64}         # labor prod own production
    p :: Array{Float64}         # hist price data
    c :: Array{Float64}         # hist cost data
    RD :: Array{Float64}        # hist R&D expenditure
    IM :: Array{Float64}        # hist immitation expenditure
    IN :: Array{Float64}        # hist innovation expenditure
    S :: Array{Float64}         # hist revenue
    HC :: Array{Float64}        # hist clients
    Π :: Array{Float64}         # hist profits
    orders :: Array             # orders
    Balance :: Balance          # balance sheet
end

"""
Checks if innovation is performed, then calls appropriate functions
"""
function innovate!(kp, global_param, all_agents, macro_struct)
    
    # determine levels of R&D, and how to divide under IN and IM
    set_RD!(kp, global_param)
    tech_choices = [(kp.A[end], kp.B[end])]

    # determine innovation of machines (Dosi et al (2010); eq. 4)
    θ_IN = 1 - exp(-global_param.ζ * kp.IN[end])
    if (θ_IN > rand())
        A_t_in = update_At(kp.A[end], global_param)
        B_t_in = update_Bt(kp.B[end], global_param)
        push!(tech_choices, (A_t_in, B_t_in))
    end

    # determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * kp.IM[end])
    if (rand(Bernoulli(θ_IM)) > rand())
        A_t_im, B_t_im = imitate_technology(kp, all_agents)
        push!(tech_choices, (A_t_im, B_t_im))
    end

    # make choice between possible technologies
    if length(tech_choices) == 1
        # if no new technologies, keep current technologies
        push!(kp.A, kp.A[end])
        push!(kp.B, kp.B[end])
    else
        c_h = map(x -> macro_struct.w[end]/x[1], tech_choices)
        p_h = map(x -> (1 + global_param.μ1)*x, c_h)
        r_h = p_h + global_param.b * c_h
        index = argmin(r_h)
        push!(kp.A, tech_choices[index][1])
        push!(kp.B, tech_choices[index][2])
        push!(kp.p, p_h[index])
        push!(kp.c, c_h[index])
    end

end


function send_brochures!(kp, all_agents, global_param)

    # set up brochure
    brochure = (kp.id, kp.p[end], kp.c[end])

    # send brochure to historical clients
    for client in kp.HC
        push!(client.brochures, brochure)
    end

    # select new clients, send brochure
    NC_potential = setdiff(all_agents.consumer_good_producers, kp.HC)

    n_choices = round(global_param.γ * length(kp.HC))
    if length(kp.HC) == 0
        n_choices = 10
    end
    
    NC = sample(NC_potential, n_choices, replace=false)
    for consumer_good_producer in NC
        push!(consumer_good_producer.brochures, brochure)
    end 

end


function imitate_technology(kp, all_agents)
    # use inverse distances as weights for choice competitor,
    distances = all_agents.capital_good_euclidian_matrix

    weights = map(x -> 1/x, distances[kp.kp_id,:])
    
    index = sample(1:length(all_agents.capital_good_producers),
                   Weights(weights))

    A_t_im = all_agents.capital_good_producers[index].A[end]
    B_t_im = all_agents.capital_good_producers[index].B[end]
    
    return A_t_im, B_t_im
end


function set_production!(kp)
    production = 0
    for order in kp.orders
        production += order[2]
    end
end


function update_At(A_last, global_param)
    # determines new labor productivity of machine produced for cp
    κ_A = rand(Beta(global_param.α1, global_param.β1))
    A_t_in = A_last*(1 + κ_A)
    return A_t_in
end


function update_Bt(B_last, global_param)
    # determines new labor productivity of own production method 
    κ_A = rand(Beta(global_param.α1, global_param.β1))
    B_t_in = B_last*(1 + κ_A)
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