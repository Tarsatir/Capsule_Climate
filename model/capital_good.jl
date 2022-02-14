using Distributions
using StatsBase


mutable struct CapitalGoodProducer <: AbstractAgent
    id :: Int                               # global id
    kp_i  :: Int                            # kp id
    A :: Vector{Float64}                     # labor prod sold product
    B :: Vector{Float64}                     # labor prod own production
    p :: Vector{Float64}                     # hist price data
    c :: Vector{Float64}                     # hist cost data
    employees :: Vector{Int}             # employees in company
    L :: Float64                            # labor units in company
    ΔLᵈ :: Float64                          # desired change in labor force
    P_FE :: Float64                         # probability of getting fired while employed
    w̄ :: Vector{Float64}                     # wage level
    wᴼ :: Float64                           # offered wage
    O :: Float64                            # total amount of machines ordered
    prod_queue :: Array                     # production queue of machines
    RD :: Vector{Float64}                    # hist R&D expenditure
    IM :: Vector{Float64}                    # hist immitation expenditure
    IN :: Vector{Float64}                    # hist innovation expenditure
    S :: Vector{Float64}                     # hist revenue
    HC :: Vector{Int}                       # hist clients
    Π :: Vector{Float64}                     # hist profits
    f :: Float64                            # market share
    brochure                                # brochure
    orders :: Array                         # orders
    Balance :: Balance                      # balance sheet
end

function initialize_kp(id :: Int, kp_i :: Int, n_captlgood :: Int, n_init_emp_kp :: Int)
    kp = CapitalGoodProducer(   # initial parameters based on rer98
        id,                     # global id
        kp_i,                   # kp id
        [1],                    # A: labor prod sold product
        [1],                    # B: labor prod own production
        [],                     # p: hist price data
        [],                     # c: hist cost data
        [],                     # employees: employees in company
        0,                      # L: labor units in company
        0,                      # ΔLᵈ: desired change in labor force
        0,                      # P_FE: probability of getting fired while employed
        [1.0],                  # w̄: wage level
        1.0,                    # wᴼ: offered wage
        0,                      # O: total amount of machines ordered
        [],                     # production queue
        [],                     # RD: hist R&D expenditure
        [],                     # IM: hist immitation expenditure
        [],                     # IN: hist innovation expenditure
        [100],                  # S: hist revenue
        [],                     # HC: hist clients
        [],                     # Π: hist profits
        1/n_captlgood,          # f: market share
        [],                     # brochure
        [],                     # orders
        Balance(               
                0.0,            # - N: inventory
                0.0,            # - K: capital
                1000.0,         # - NW: liquid assets
                0.0,            # - Deb: debt
                0.0             # - EQ: equity
            )               
    )
    return kp
end


"""
Checks if innovation is performed, then calls appropriate functions
"""
function innovate_kp!(
    kp::CapitalGoodProducer, 
    global_param, 
    all_kp::Vector{Int}, 
    kp_distance_matrix::Array{Float64}, 
    model::ABM
    )
    
    # determine levels of R&D, and how to divide under IN and IM
    set_RD_kp!(kp, global_param.ξ, global_param.ν)
    tech_choices = [(kp.A[end], kp.B[end])]

    # determine innovation of machines (Dosi et al (2010); eq. 4)
    θ_IN = 1 - exp(-global_param.ζ * kp.IN[end])
    if (rand(Bernoulli(θ_IN)) > rand())
        A_t_in = update_At_kp(kp.A[end], global_param)
        B_t_in = update_Bt_kp(kp.B[end], global_param)
        push!(tech_choices, (A_t_in, B_t_in))
    end

    # determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * kp.IM[end])
    if (rand(Bernoulli(θ_IM)) > rand())
        A_t_im, B_t_im = imitate_technology_kp(kp, all_kp, kp_distance_matrix, model)
        push!(tech_choices, (A_t_im, B_t_im))
    end

    # make choice between possible technologies
    if length(tech_choices) == 1
        # if no new technologies, keep current technologies
        push!(kp.A, kp.A[end])
        push!(kp.B, kp.B[end])
        c_h = kp.w̄[end]/kp.B[end]
        p_h = (1 + global_param.μ1) * c_h
        push!(kp.c, c_h)
        push!(kp.p, p_h)
    else
        # if new technologies, update price data
        c_h = map(x -> kp.w̄[end]/x[1], tech_choices)
        p_h = map(x -> (1 + global_param.μ1)*x, c_h)
        r_h = p_h + global_param.b * c_h
        index = argmin(r_h)
        push!(kp.A, tech_choices[index][1])
        push!(kp.B, tech_choices[index][2])
        push!(kp.c, c_h[index])
        push!(kp.p, p_h[index])
    end

end


function send_brochures_kp!(
    kp::CapitalGoodProducer,
    all_cp::Vector{Int}, 
    global_param,
    model::ABM
    )

    # set up brochure
    brochure = (kp.id, kp.p[end], kp.c[end], kp.A[end])
    kp.brochure = brochure

    # send brochure to historical clients
    for cp_id in kp.HC
        cp = model[cp_id]
        push!(cp.brochures, brochure)
    end

    # select new clients, send brochure
    NC_potential = setdiff(all_cp, kp.HC)

    n_choices = Int(round(global_param.γ * length(kp.HC)))
    
    # send brochures to new clients
    NC = sample(NC_potential, n_choices; replace=false)
    for cp_id in NC
        cp = model[cp_id]
        push!(cp.brochures, brochure)
    end 

end


"""

Uses inverse distances as weights for choice competitor to immitate
"""
function imitate_technology_kp(
    kp::CapitalGoodProducer, 
    all_kp::Vector{Int}, 
    kp_distance_matrix, 
    model::ABM
    ) :: Tuple{Float64, Float64}

    weights = map(x -> 1/x, kp_distance_matrix[kp.kp_i,:])
    
    idx = sample(all_kp, Weights(weights))

    A_t_im = model[idx].A[end]
    B_t_im = model[idx].B[end]
    
    return A_t_im, B_t_im
end


function update_At_kp(
    A_last :: Float64, 
    global_param
    )
    # determines new labor productivity of machine produced for cp
    κ_A = min(rand(Beta(global_param.α1, global_param.β1)), global_param.κ_upper)
    A_t_in = A_last*(1 + κ_A)
    return A_t_in
end


function update_Bt_kp(
    B_last :: Float64, 
    global_param
    )
    # determines new labor productivity of own production method 
    κ_B = min(rand(Beta(global_param.α1, global_param.β1)), global_param.κ_upper)
    # println("κ ", κ_B)
    B_t_in = B_last*(1 + κ_B)
    return B_t_in
end


function set_price_kp!(
    kp::CapitalGoodProducer, 
    μ1::Float64, 
    )
    c_t = kp.w̄[end] / kp.B[end]
    p_t = (1 + μ1) * c_t
    push!(kp.c, c_t)
    push!(kp.p, p_t)
end


"""
Determines the level of R&D, and how it is divided under innovation (IN) 
and immitation (IM). based on Dosi et al (2010)
"""
function set_RD_kp!(
    kp::CapitalGoodProducer, 
    ξ::Float64, 
    ν::Float64
    )
    # determine total R&D spending at time t, (Dosi et al, 2010; eq. 3)
    RD_new = ν * kp.S[end]
    push!(kp.RD, RD_new)

    # decide fractions innovation (IN) and immitation (IM), 
    # (Dosi et al, 2010; eq. 3.5)
    IN_t_new = ξ * RD_new
    IM_t_new = (1 - ξ) * RD_new
    push!(kp.IN, IN_t_new)
    push!(kp.IM, IM_t_new)
end


function plan_production_kp!(
    kp::CapitalGoodProducer
    )
    
    if (length(kp.orders) > 0)
        # determine total amount of machines to produce
        kp.O = sum(map(order -> order[2], kp.orders))
    end

    # determine amount of labor to hire
    kp.ΔLᵈ = kp.O/kp.B[end] - kp.L
end


function produce_goods_kp!(
    kp::CapitalGoodProducer
    )

    # println(kp.B[end] * kp.L, " ", kp.O)

    # check if enough labor available for all machines
    unocc_L = kp.L
    for order in kp.orders
        req_L = order[2] / kp.B[end]
        if unocc_L > req_L
            push!(kp.prod_queue, order)
            unocc_L -= req_L
        end
    end

    # reset order amount back to zero
    kp.O = 0
end


"""
Sends orders from production queue to cp
"""
function send_orders_kp!(
    kp::CapitalGoodProducer,
    model::ABM
    )

    tot_freq = 0

    for order in kp.prod_queue

        cp_id = order[1]
        Iₜ = order[2]

        machine = initialize_machine()
        machine.A = kp.A[end]
        machine.freq = Iₜ
        # TODO check if this is correct

        tot_freq += machine.freq

        receive_machines!(model[cp_id], machine, Iₜ)
    end
    
    S = tot_freq * kp.p[end]
    Π = tot_freq * (kp.p[end] - kp.c[end])

    # println(Π, " ", tot_freq)

    push!(kp.S, S)
    push!(kp.Π, Π)

    # TODO: request money for machines that were produced

    # TODO: describe how labor market frictions affect delivery machines

    # Empty production queue
    reset_order_queue_kp!(kp)

end


function reset_order_queue_kp!(
    kp:: CapitalGoodProducer
    )
    kp.orders = []
end


function select_HC_kp!(
    kp::CapitalGoodProducer, 
    all_cp::Vector{Int}
    )
    kp.HC = sample(all_cp, 10; replace=false)
end