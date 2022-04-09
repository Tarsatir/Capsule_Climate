@with_kw mutable struct CapitalGoodProducer <: AbstractAgent
    # T::Int
    id::Int                               # global id
    kp_i::Int                             # kp id, used for distance matrix
    age::Int = 0                          # firm age
    first_gen::Bool                       # shows if producer is in first generation
    
    # Technology and innovation
    A::Float64 = 1.0                      # labor prod sold product
    B::Float64 = 1.0                      # labor prod own production
    RD::Float64 = 0.0                     # R&D expenditure
    IM::Float64 = 0.0                     # hist immitation expenditure
    IN::Float64 = 0.0                     # hist innovation expenditure

    # Price and cost data
    μ::Vector{Float64}                    # markup rates
    p::Vector{Float64} = fill(1+μ[end], 3) # hist price data
    c::Vector{Float64} = ones(Float64, 3) # hist cost data
    
    # Employment
    employees::Vector{Int} = []           # employees in company
    L::Float64 = 0.0                      # labor units in company
    ΔLᵈ::Float64 = 0.0                    # desired change in labor force
    w̄::Vector{Float64}                    # wage level
    wᴼ::Float64 = w̄[end]                  # offered wage
    wᴼ_max::Float64 = 0.0                 # maximum offered wage

    O::Float64 = 0.0                      # total amount of machines ordered
    prod_queue::Array = []                # production queue of machines
    Q::Vector{Float64} = zeros(Float64, 3)# production quantities
    D::Vector{Float64} = zeros(Float64, 3)# hist demand
    HC::Vector{Int} = []                  # hist clients
    Π::Vector{Float64} = zeros(Float64, 3)# hist profits
    Πᵀ::Vector{Float64} = zeros(Float64, 3)# hist profits after tax
    debt_installments::Vector{Float64} = zeros(Float64, 4)   # installments of debt repayments
    f::Vector{Float64}                    # market share
    brochure = []                         # brochure
    orders::Array = []                    # orders
    balance::Balance                      # balance sheet
    curracc::FirmCurrentAccount = FirmCurrentAccount() # current account
end


"""
Initializes kp agent, default is the heterogeneous state, otherwise properties are given
    as optional arguments.
"""
function initialize_kp(
    id::Int, 
    kp_i::Int,
    n_captlgood::Int;
    NW=1000,
    A=1.0,
    B=1.0,
    p=1.0,
    μ=0.2,
    w̄=1.0,
    wᴼ=1.0,
    Q=100,
    D=100,
    f=1/n_captlgood,
    first_gen=true
    )

    # balance = Balance(               
    #     0.0,                    # - N: inventory
    #     0.0,                    # - K: capital
    #     NW,                     # - NW: liquid assets
    #     0.0,                    # - debt: debt
    #     NW                      # - EQ: equity
    # )
    
    # curracc = FirmCurrentAccount(0,0,0,0,0,0,0)

    # kp = CapitalGoodProducer(   
    #         id = id,                     # global id
    #         kp_i,                   # kp_i, used for distance matrix
    #         # 0,                      # firm age
    #         first_gen,              # shows if producer is in first generation
    #         [A],                    # A: labor prod sold product
    #         [B],                    # B: labor prod own production
    #         [p],                    # p: hist price data
    #         [],                     # c: hist cost data
    #         [μ],                    # μ: hist markup rate
    #         [],                     # employees: employees in company
    #         0,                      # L: labor units in company
    #         0,                      # ΔLᵈ: desired change in labor force
    #         [w̄],                    # w̄: wage level
    #         wᴼ,                     # wᴼ: offered wage
    #         0.0,                    # wᴼ_max: expected offered wage
    #         0,                      # O: total amount of machines ordered
    #         [],                     # production queue
    #         [],                     # Q: production quantity
    #         [],                     # RD: hist R&D expenditure
    #         [],                     # IM: hist immitation expenditure
    #         [],                     # IN: hist innovation expenditure
    #         [],                     # D: hist demand
    #         [],                     # HC: hist clients
    #         [0.0],                  # Π: hist profits
    #         zeros(4),               # debt installments
    #         [f],                    # f: market share
    #         [],                     # brochure
    #         [],                     # orders
    #         balance,                # balance
    #         curracc                 # current account   
    #     )

        kp = CapitalGoodProducer(
            id = id,                        
            kp_i = kp_i,                   
            # 0,                      
            first_gen = first_gen,          

            A = A,                        
            B = B,                       
            # RD = zeros(Float64, T),
            # IM = zeros(Float64, T),
            # IN = zeros(Float64, T),

            μ = fill(μ, 3),
            # p = fill(p, T),                 
            # c = ones(Float64, T),

            # [],                     # employees: employees in company
            # 0,                      # L: labor units in company
            # 0,                      # ΔLᵈ: desired change in labor force
            w̄ = fill(w̄, 3), 
            wᴼ = w̄,                     # wᴼ: offered wage
            # 0.0,                    # wᴼ_max: expected offered wage

            # 0,                      # O: total amount of machines ordered
            # [],                     # production queue
            # Q = zeros(Float64, T),   
            # [],                     # RD: hist R&D expenditure
            # [],                     # IM: hist immitation expenditure
            # [],                     # IN: hist innovation expenditure
            # D = zeros(Float64, T),  
            # [],                     # HC: hist clients
            # Π = zeros(Float64, T),    # Π: hist profits
            # Πᵀ = zeros(Float64, T),
            # zeros(4),               # debt installments
            f = fill(f, 3),           # f: market share
            # [],                     # brochure
            # [],                     # orders
            balance = Balance(NW=NW, EQ=NW),                # balance
            # curracc                 # current account   
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
    tech_choices = [(kp.A, kp.B)]

    # Determine innovation of machines (Dosi et al (2010); eq. 4)
    θ_IN = 1 - exp(-global_param.ζ * kp.IN)
    if rand(Bernoulli(θ_IN))
        A_t_in = update_At_kp(kp.A, global_param)
        B_t_in = update_Bt_kp(kp.B, global_param)
        if A_t_in > kp.A || B_t_in > kp.B
            push!(tech_choices, (A_t_in, B_t_in))
        end
    end

    # TODO compute real value innovation like rer98

    # Determine immitation of competitors
    θ_IM = 1 - exp(-global_param.ζ * kp.IM)
    if rand(Bernoulli(θ_IM))
        A_t_im, B_t_im = imitate_technology_kp(kp, all_kp, kp_distance_matrix, model)
        if A_t_im > kp.A || B_t_im > kp.B
            push!(tech_choices, (A_t_im, B_t_im))
        end
    end

    choose_technology_kp!(kp, w̄, global_param, tech_choices)
end


"""
Lets kp choose technology
"""
function choose_technology_kp!(
    kp::CapitalGoodProducer,
    w̄::Float64,
    global_param::GlobalParam,
    tech_choices
    )

    # TODO: DESCRIBE
    update_μ_kp!(kp)

    # Make choice between possible technologies
    if length(tech_choices) == 1
        # If no new technologies, keep current technologies
        # push!(kp.A, kp.A[end])
        # kp.A[t] = kp.A[t-1]
        # push!(kp.B, kp.B[end])
        # kp.B[t] = kp.B[t-1]

        c_h = kp.w̄[end]/kp.B
        p_h = (1 + kp.μ[end]) * c_h
        # push!(kp.c, c_h)
        # push!(kp.p, p_h)
        shift_and_append!(kp.c, c_h)
        shift_and_append!(kp.p, p_h)
    else
        # If new technologies, update price data
        c_h_cp = map(tech -> w̄/tech[1], tech_choices)
        c_h_kp = map(tech -> kp.w̄[end]/tech[2], tech_choices)
 
        p_h = map(c -> (1 + kp.μ[end])*c, c_h_kp)
        r_h = p_h + global_param.b * c_h_cp 
        idx = argmin(r_h)

        kp.A = tech_choices[idx][1]
        kp.B = tech_choices[idx][2]

        # push!(kp.A, tech_choices[idx][1])
        # push!(kp.B, tech_choices[idx][2])
        # push!(kp.c, c_h_kp[idx])
        # push!(kp.p, p_h[idx])
        shift_and_append!(kp.c, c_h_kp[idx])
        shift_and_append!(kp.p, p_h[idx])
    end
end


"""
Creates brochures and sends to potential clients.
"""
function send_brochures_kp!(
    kp::CapitalGoodProducer,
    all_cp::Vector{Int}, 
    global_param,
    model::ABM;
    n_hist_clients=50::Int
    )

    # Set up brochure
    brochure = (kp.id, kp.p[end], kp.c[end], kp.A)
    kp.brochure = brochure

    # Send brochure to historical clients
    for cp_id in kp.HC
        push!(model[cp_id].brochures, brochure)
    end

    # Select new clients, send brochure
    NC_potential = setdiff(all_cp, kp.HC)

    if length(kp.HC) == 0
        n_choices = n_hist_clients
    else
        n_choices = Int(round(global_param.γ * length(kp.HC)))
    end
    
    # Send brochures to new clients
    NC = sample(NC_potential, min(n_choices, length(NC_potential)); replace=false)
    for cp_id in NC
        push!(model[cp_id].brochures, brochure)
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

    weights = map(x -> 1/x, kp_distance_matrix[kp.kp_i,:])
    idx = sample(all_kp, Weights(weights))

    A_t_im = model[idx].A
    B_t_im = model[idx].B
    
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

    c_t = kp.w̄[end] / kp.B
    # p_t = (1 + μ1) * c_t
    p_t = (1 + kp.μ[end]) * c_t
    # push!(kp.c, c_t)
    # push!(kp.p, p_t)
    shift_and_append!(kp.c, c_t)
    shift_and_append!(kp.p, p_t)
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
    # TODO: describe change of added NW
    
    if length(kp.Q) > 0
        prev_S = kp.p[end] * kp.Q[end]
        RD_new = ν * prev_S
        if prev_S == 0
            RD_new = 0.3 * max(kp.balance.NW, 0)
        end
    else
        RD_new = 0.3 * max(kp.balance.NW, 0)
    end

    # TODO: now based on prev profit to avoid large losses. If kept, describe!

    # push!(kp.RD, RD_new)
    kp.RD = RD_new
    kp.curracc.TCI += RD_new

    # Decide fractions innovation (IN) and immitation (IM), 
    # (Dosi et al, 2010; eq. 3.5)
    kp.IN = ξ * RD_new
    kp.IM = (1 - ξ) * RD_new
    # push!(kp.IN, IN_t_new)
    # push!(kp.IM, IM_t_new)
end


"""
Lets kp receive orders, adds client as historical clients if it is not yet.
"""
function receive_order_kp!(
    kp::CapitalGoodProducer,
    cp_id::Int
    )

    push!(kp.orders, cp_id)

    # If cp not in HC yet, add as a historical client
    if cp_id ∉ kp.HC
        push!(kp.HC, cp_id)
    end
end


"""
Based on received orders, sets labor demand to fulfill production.
"""
function plan_production_kp!(
    kp::CapitalGoodProducer,
    global_param::GlobalParam,
    model::ABM
    )

    update_w̄_p!(kp, model)
    
    # Determine total amount of capital units to produce and amount of labor to hire
    kp.O = length(kp.orders) * global_param.freq_per_machine

    # Determine amount of labor to hire
    kp.ΔLᵈ = kp.O / kp.B + kp.RD / kp.w̄[end] - kp.L

    update_wᴼ_max_kp!(kp)
end


"""
Lets kp add goods to the production queue, based on available labor supply
"""
function produce_goods_kp!(
    kp::CapitalGoodProducer,
    global_param::GlobalParam,
    )

    # Determine what the total demand is, regardless if it can be satisfied
    D = length(kp.orders) * global_param.freq_per_machine
    push!(kp.D, D)

    # Determine how much labor is needed to produce a full machine
    req_L = global_param.freq_per_machine / kp.B

    # Check if production is constrained
    if kp.L >= req_L * length(kp.orders)

        # Enough labor available, perform full production
        kp.prod_queue = kp.orders

    else

        # Production constrained, determine amount of production possible
        # and randomly select which machines to produce
        n_poss_prod = floor(Int, kp.L / req_L)
        kp.prod_queue = sample(kp.orders, n_poss_prod; replace=false)

    end

    # Append total production amount of capital units
    Q = length(kp.prod_queue) * global_param.freq_per_machine
    push!(kp.Q, Q)

    # println("len orders: $(length(kp.orders)), len prod queue: $(length(kp.prod_queue))")

    # Empty order queue
    kp.orders = []
end


"""
Sends orders from production queue to cp.
"""
function send_orders_kp!(
    kp::CapitalGoodProducer,
    global_param::GlobalParam,
    model::ABM
    )

    if length(kp.prod_queue) == 0
        return nothing
    end

    # Count how many machines each individual cp ordered
    # println(kp.prod_queue)
    machines_per_cp = counter(kp.prod_queue)
    # println(machines_per_cp)

    for (cp_id, n_machines) in machines_per_cp

        # Produce machines in production queue, send to cp
        machines = initialize_machine_stock(global_param.freq_per_machine, n_machines;
                                            p=kp.p[end], A=kp.A)
        Iₜ = n_machines * global_param.freq_per_machine * kp.p[end]
        receive_machines_cp!(model[cp_id], machines, Iₜ)
    end
    
    # Update sales
    kp.curracc.S = kp.Q[end] * kp.p[end]

    # Empty production queue
    kp.prod_queue = []
end


function reset_order_queue_kp!(
    kp::CapitalGoodProducer
    )

    kp.orders = []
end


function select_HC_kp!(
    kp::CapitalGoodProducer, 
    all_cp::Vector{Int};
    n_hist_clients=10::Int
    )

    kp.HC = sample(all_cp, n_hist_clients; replace=false)
end


function update_wᴼ_max_kp!(
    kp::CapitalGoodProducer
    )
    # TODO: DESCRIBE IN MODEL
    kp.wᴼ_max = kp.B * kp.p[end] 
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


"""
Updates market share of all kp.
"""
function update_marketshare_kp!(
    all_kp::Vector{Int},
    model::ABM
    )

    kp_market = sum(kp_id -> model[kp_id].D[end], all_kp)

    for kp_id in all_kp
        if kp_market == 0
            f = 1 / length(all_kp)
        else
            f = model[kp_id].D[end] / kp_market
        end
        push!(model[kp_id].f, f)
    end
end


# TRIAL: DESCRIBE
function update_μ_kp!(
    kp::CapitalGoodProducer
    )

    # b = 0.3
    # l = 2

    # if length(kp.Π) > l
    #     # mean_μ = mean(kp.μ[end-l:end-1])
    #     # Δμ = (cp.μ[end] - cp.μ[end-1]) / cp.μ[end-1]
    #     # Δμ = (kp.μ[end] - mean_μ) / mean_μ

    #     # mean_Π = mean(kp.Π[end-l:end-1])
    #     # ΔΠ = (cp.Π[end] - cp.Π[end-1]) / cp.Π[end-1]
    #     # ΔΠ = (kp.Π[end] - mean_Π) / mean_Π
    #     shock_sign = 1
    #     if kp.Π[end] <= kp.Π[end-1]
    #         shock_sign = -1
    #     end

    #     # println("$mean_μ, $mean_Π")
    #     # println("Δμ: $Δμ, $(sign(Δμ)), ΔΠ: $ΔΠ, $(sign(ΔΠ))")

    #     shock = rand(Uniform(0.0, b))

    #     new_μ = max(kp.μ[end] * (1 + shock_sign * shock), 0)
    #     push!(kp.μ, new_μ)

    # elseif kp.Π[end] == 0
    #     push!(kp.μ, kp.μ[end] * (1 + rand(Uniform(-b, 0.0))))
    # else
    #     push!(kp.μ, kp.μ[end] * (1 + rand(Uniform(-b, b))))
    # end

    push!(kp.μ, kp.μ[end])

end


"""
Replaces bankrupt kp with new kp. Gives them a level of technology and expectations
    from another kp. 
"""
function replace_bankrupt_kp!(
    bankrupt_kp::Vector{Int},
    bankrupt_kp_i::Vector{Int},
    all_kp::Vector{Int},
    global_param::GlobalParam,
    indexfund_struct::IndexFund,
    init_param::InitParam,
    t::Int,
    model::ABM
    )

    # TODO: describe in model

    # Check if not all kp have gone bankrupt, in this case, 
    # kp with highest NW will not be removed
    if length(bankrupt_kp) == length(all_kp)
        kp_id_max_NW = all_kp[1]
        for i in 2:length(all_kp)
            if model[all_kp[i]].curracc.NW > model[kp_id_max_NW].curracc.NW
                kp_id_max_NW = all_kp[i]
            end
        end
        bankrupt_kp = bankrupt_kp[!kp_id_max_NW]
    end

    # Determine all possible kp and set weights for sampling proportional to the 
    # quality of their technology
    poss_kp = filter(kp_id -> kp_id ∉ bankrupt_kp, all_kp)
    # weights = map(kp_id -> min(model[kp_id].A, model[kp_id].B), poss_kp)

    # Get the technology frontier
    A_max = 0
    A_min = 0
    A_max_id = 0

    B_max = 0
    B_min = 0

    for kp_id in poss_kp
        # Check if A is max or min in population
        if model[kp_id].A > A_max
            A_max = model[kp_id].A
            A_max_id = kp_id
        elseif model[kp_id].A < A_min
            A_min = model[kp_id].A
        end

        # Check if B is max or min in population
        if model[kp_id].B > B_max
            B_max = model[kp_id].B
        elseif model[kp_id].B < B_min
            B_min = model[kp_id].B_min
        end 
    end

    # Compute the average stock of liquid assets of non-bankrupt kp
    avg_NW = mean(kp_id -> model[kp_id].balance.NW, poss_kp)
    NW_coefficients = rand(Uniform(global_param.φ3, global_param.φ4),
                           length(bankrupt_kp))

    # Compute share of investments that can be paid from the investment fund                       
    all_req_NW = sum(avg_NW .* NW_coefficients)
    frac_NW_if = decide_investments_if!(indexfund_struct, all_req_NW, t)

    # Re-use id of bankrupted company
    for (i, (kp_id, kp_i)) in enumerate(zip(bankrupt_kp, bankrupt_kp_i))
        # Sample a producer of which to take over the technologies, proportional to the 
        # quality of the technology
        tech_coeff = (global_param.φ5 + rand(Beta(global_param.α2, global_param.β2)) 
                                        * (global_param.φ6 - global_param.φ5))
        new_A = max(A_max * (1 + tech_coeff), A_min, init_param.A_0)
        new_B = max(B_max * (1 + tech_coeff), B_min, init_param.B_0)

        NW_stock = NW_coefficients[i] * avg_NW

        # Initialize new kp
        new_kp = initialize_kp(
            kp_id, 
            kp_i, 
            length(all_kp);
            NW=NW_stock,
            A=new_A,
            B=new_B,
            μ=model[A_max_id].μ[end],
            w̄=model[A_max_id].w̄[end],
            f=0.0,
            first_gen=false
        )

        # Borrow the remaining funds
        borrow_funds_p!(new_kp, (1 - frac_NW_if) * NW_stock, global_param.b)

        add_agent!(new_kp, model)
    end
end