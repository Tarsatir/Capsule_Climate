mutable struct CapitalGoodProducer <: AbstractAgent
    id :: Int                               # global id
    kp_i :: Int                             # kp id, used for distance matrix
    A :: Vector{Float64}                    # labor prod sold product
    B :: Vector{Float64}                    # labor prod own production
    p :: Vector{Float64}                    # hist price data
    c :: Vector{Float64}                    # hist cost data
    employees :: Vector{Int}                # employees in company
    L :: Float64                            # labor units in company
    ΔLᵈ :: Float64                          # desired change in labor force
    P_FE :: Float64                         # probability of getting fired while employed
    w̄ :: Vector{Float64}                    # wage level
    wᴼ :: Float64                           # offered wage
    wᴼ_max :: Float64                       # maximum offered wage
    O :: Float64                            # total amount of machines ordered
    prod_queue :: Array                     # production queue of machines
    Q :: Vector{Float64}                    # production quantities
    RD :: Vector{Float64}                   # hist R&D expenditure
    IM :: Vector{Float64}                   # hist immitation expenditure
    IN :: Vector{Float64}                   # hist innovation expenditure
    D :: Vector{Float64}                    # hist revenue
    HC :: Vector{Int}                       # hist clients
    Π :: Vector{Float64}                    # hist profits
    debt_installments :: Vector{Float64}   # installments of debt repayments
    f :: Float64                            # market share
    brochure                                # brochure
    orders :: Array                         # orders
    balance :: Balance                      # balance sheet
    curracc :: FirmCurrentAccount           # current account
end

function initialize_kp(
    id::Int, 
    n_captlgood::Int,
    kp_i::Int,
    n_init_emp_kp::Int
    )

    balance = Balance(               
        0.0,                    # - N: inventory
        0.0,                    # - K: capital
        1000.0,                 # - NW: liquid assets
        0.0,                    # - debt: debt
        0.0                     # - EQ: equity
    )
    
    curracc = FirmCurrentAccount(0,0,0,0,0,0,0)

    kp = CapitalGoodProducer(   # initial parameters based on rer98
        id,                     # global id
        kp_i,                   # kp_i, used for distance matrix
        [1],                    # A: labor prod sold product
        [1],                    # B: labor prod own production
        [1.0],                  # p: hist price data
        [],                     # c: hist cost data
        [],                     # employees: employees in company
        0,                      # L: labor units in company
        0,                      # ΔLᵈ: desired change in labor force
        0,                      # P_FE: probability of getting fired while employed
        [1.0],                  # w̄: wage level
        1.0,                    # wᴼ: offered wage
        1.0,                    # wᴼₑ: expected offered wage
        0,                      # O: total amount of machines ordered
        [],                     # production queue
        [10],                   # Q: production quantity
        [],                     # RD: hist R&D expenditure
        [],                     # IM: hist immitation expenditure
        [],                     # IN: hist innovation expenditure
        [100],                  # D: hist revenue
        [],                     # HC: hist clients
        [0.0],                  # Π: hist profits
        fill(0.0, 4),           # debt installments
        1/n_captlgood,          # f: market share
        [],                     # brochure
        [],                     # orders
        balance,                # balance
        curracc                 # current account   
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
    w̄::Float64,
    model::ABM,
    )
    
    # TODO include labor demand for researchers

    # Determine levels of R&D, and how to divide under IN and IM
    set_RD_kp!(kp, global_param.ξ, global_param.ν)
    tech_choices = [(kp.A[end], kp.B[end])]

    # Determine innovation of machines (Dosi et al (2010); eq. 4)
    θ_IN = 1 - exp(-global_param.ζ * kp.IN[end])
    if rand(Bernoulli(θ_IN))
        A_t_in = update_At_kp(kp.A[end], global_param)
        B_t_in = update_Bt_kp(kp.B[end], global_param)
        if A_t_in > kp.A[end] || B_t_in > kp.B[end]
            push!(tech_choices, (A_t_in, B_t_in))
        end
    end

    # TODO compute real value innovation like rer98

    # Determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * kp.IM[end])
    if rand(Bernoulli(θ_IM))
        A_t_im, B_t_im = imitate_technology_kp(kp, all_kp, kp_distance_matrix, model)
        if A_t_im > kp.A[end] || B_t_im > kp.B[end]
            push!(tech_choices, (A_t_im, B_t_im))
        end
    end

    choose_technology_kp!(kp, w̄, tech_choices)
end


"""
Lets kp choose technology
"""
function choose_technology_kp!(
    kp::CapitalGoodProducer,
    w̄::Float64,
    tech_choices
    )

    # Make choice between possible technologies
    if length(tech_choices) == 1
        # If no new technologies, keep current technologies
        push!(kp.A, kp.A[end])
        push!(kp.B, kp.B[end])
        c_h = kp.w̄[end]/kp.B[end]
        p_h = (1 + global_param.μ1) * c_h
        push!(kp.c, c_h)
        push!(kp.p, p_h)
    else
        # If new technologies, update price data
        c_h_kp = map(tech -> kp.w̄[end]/tech[2], tech_choices)
        c_h_cp = map(tech -> w̄/tech[1], tech_choices)
        p_h = map(c -> (1 + global_param.μ1)*c, c_h_kp)
        r_h = p_h + global_param.b * c_h_cp 
        index = argmin(r_h)
        push!(kp.A, tech_choices[index][1])
        push!(kp.B, tech_choices[index][2])
        push!(kp.c, c_h_kp[index])
        push!(kp.p, p_h[index])
    end
end


"""
Creates brochures and sends to potential clients.
"""
function send_brochures_kp!(
    kp::CapitalGoodProducer,
    all_cp::Vector{Int}, 
    global_param,
    model::ABM
    )

    # Set up brochure
    brochure = (kp.id, kp.p[end], kp.c[end], kp.A[end])
    kp.brochure = brochure

    # Send brochure to historical clients
    for cp_id in kp.HC
        cp = model[cp_id]
        push!(cp.brochures, brochure)
    end

    # Select new clients, send brochure
    NC_potential = setdiff(all_cp, kp.HC)

    n_choices = Int(round(global_param.γ * length(kp.HC)))
    
    # Send brochures to new clients
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
    )::Tuple{Float64, Float64}

    # weights = map(x -> 1/x, kp_distance_matrix[kp.kp_i,:])
    # idx = sample(all_kp, Weights(weights))
    # TODO: uncomment

    idx = sample(all_kp)

    A_t_im = model[idx].A[end]
    B_t_im = model[idx].B[end]
    
    return A_t_im, B_t_im
end


"""
Determines new labor productivity of machine produced for cp
"""
function update_At_kp(
    A_last::Float64, 
    global_param::GlobalParam
    )::Float64
    
    κ_A = rand(Beta(global_param.α1, global_param.β1))
    κ_A = global_param.κ_lower + κ_A * (global_param.κ_upper - global_param.κ_lower)

    A_t_in = A_last * (1 + κ_A)
    return A_t_in
end


"""
Determines new labor productivity of own production method.
"""
function update_Bt_kp(
    B_last::Float64, 
    global_param::GlobalParam
    )::Float64

    κ_B = rand(Beta(global_param.α1, global_param.β1))
    κ_B = global_param.κ_lower + κ_B * (global_param.κ_upper - global_param.κ_lower)

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

    # Determine total R&D spending at time t, (Dosi et al, 2010; eq. 3)
    prev_S = kp.p[end] * kp.Q[end]
    RD_new = ν * prev_S
    # RD_new = max(ν * kp.Π[end], 0)

    # TODO: now based on prev profit to avoid large losses. If kept, describe!

    push!(kp.RD, RD_new)
    kp.curracc.TCI += RD_new

    # Decide fractions innovation (IN) and immitation (IM), 
    # (Dosi et al, 2010; eq. 3.5)
    IN_t_new = ξ * RD_new
    IM_t_new = (1 - ξ) * RD_new
    push!(kp.IN, IN_t_new)
    push!(kp.IM, IM_t_new)
end


"""
Based on received orders, sets labor demand to fulfill production.
"""
function plan_production_kp!(
    kp::CapitalGoodProducer,
    model::ABM
    )

    update_w̄_p!(kp, model)
    
    # Determine total amount of machines to produce
    if length(kp.orders) > 0
        kp.O = sum(map(order -> order[2] / kp.p[end], kp.orders)) 
    else
        kp.O = 0.0
    end

    # Determine amount of labor to hire
    kp.ΔLᵈ = 50 * floor(max(kp.O / kp.B[end] + kp.RD[end] / kp.w̄[end] - kp.L, -kp.L) / 50)
    # kp.ΔLᵈ = max(kp.O / kp.B[end] + kp.RD[end] / kp.w̄[end] - kp.L, -kp.L)

    # println("labor ", kp.ΔLᵈ, " ", kp.O / kp.B[end], " ", kp.RD[end] / kp.w̄[end], " ", kp.L)
    # kp.ΔLᵈ = kp.O / kp.B[end] - kp.L

    update_wᴼ_max_kp!(kp)
end


"""
Lets kp produce ordered goods based on actual available labor stock.
"""
function produce_goods_kp!(
    kp::CapitalGoodProducer
    )

    total_Q = 0

    # Check if enough labor available for all machines
    unocc_L = kp.L
    for order in kp.orders

        # Only push order in production queue if enough labor available
        q = order[2] / kp.p[end]
        req_L = q / kp.B[end]

        if unocc_L > req_L
            # Full order can be fulfilled
            total_Q += q
            push!(kp.prod_queue, order)
            unocc_L -= req_L
        elseif unocc_L > 0.0
            # Partial order can be fulfilled
            part_q = unocc_L * kp.B[end]
            total_Q += part_q
            new_Iₜ = part_q * kp.p[end]
            push!(kp.prod_queue, (order[1], new_Iₜ))
            unocc_L = 0.0
        end
    end

    # Reset order amount back to zero
    # kp.O = 0

    # println(length(kp.orders), " ", length(kp.prod_queue), " ", kp.L, " ", kp.O)

    # Push production amount
    push!(kp.Q, total_Q)
end


"""
Sends orders from production queue to cp.
"""
function send_orders_kp!(
    kp::CapitalGoodProducer,
    model::ABM
    )

    tot_freq = 0
    for order in kp.prod_queue

        cp_id = order[1]
        Iₜ = order[2]
        freq = Iₜ / kp.p[end]

        machine = initialize_machine(freq, 0, kp.p[end], kp.A[end])

        tot_freq += freq

        receive_machines!(model[cp_id], machine, Iₜ)
    end
    
    # Add demand and sales.
    # D = tot_freq
    # push!(kp.D, D)

    kp.curracc.S = tot_freq * kp.p[end]

    # TODO: request money for machines that were produced

    # TODO: describe how labor market frictions affect delivery machines

    # Empty order and production queue
    reset_order_queue_kp!(kp)
    reset_prod_queue_kp!(kp)
end


function reset_order_queue_kp!(
    kp::CapitalGoodProducer
    )

    kp.orders = []
end


function reset_prod_queue_kp!(
    kp::CapitalGoodProducer
    )

    kp.prod_queue = []
end


function select_HC_kp!(
    kp::CapitalGoodProducer, 
    all_cp::Vector{Int}
    )

    kp.HC = sample(all_cp, 10; replace=false)
end


function update_wᴼ_max_kp!(
    kp::CapitalGoodProducer
    )
    # TODO: DESCRIBE IN MODEL
    kp.wᴼ_max = kp.B[end] * kp.p[end] 
    # if kp.ΔLᵈ > 0
    #     kp.wᴼ_max = (kp.O * kp.p[end] - kp.w̄[end] * kp.L) / kp.ΔLᵈ
    # else
    #     kp.wᴼ_max = 0
    # end
end


"""
Filters out historical clients if they went bankrupt
"""
function remove_bankrupt_HC_kp!(
    kp::CapitalGoodProducer,
    bankrupt_lp::Vector{Int},
    bankrupt_bp::Vector{Int}
    )

    filter!(bp_id -> bp_id ∉ bankrupt_bp, kp.HC)
    filter!(lp_id -> lp_id ∉ bankrupt_lp, kp.HC)
end