using Statistics


"""
Defines struct for consumer good producer
"""
mutable struct ConsumerGoodProducer <: AbstractAgent
    id :: Int                   # id
    type_good :: String         # type of good produced by producer 
    p :: Vector{Float64}        # hist prices
    c :: Vector{Float64}        # hist cost
    RD :: Vector{Float64}       # R&D spending
    D :: Vector{Float64}        # hist demand
    Dᵉ :: Float64               # exp demand
    hh_queue :: Vector          # vector containing orders of demanding households
    Nᵈ :: Float64               # desired inventory
    Q :: Vector{Float64}        # hist production
    Qᵉ :: Float64               # exp production
    I :: Vector{Float64}        # hist investments
    Ξ :: Vector{Machine}        # capital stock
    K :: Float64                # total amount of machines
    RS :: Vector{Machine}       # list of to-be replaced machines
    employees:: Vector{Int}     # employees list
    L :: Float64                # labor units
    Lᵉ:: Float64                # exp labor force
    ΔLᵈ :: Float64              # desired change in labor force
    # P_FE :: Float64             # probability of getting fired while employed
    w̄ :: Vector{Float64}        # wage level
    wᴼ :: Float64               # offered wage
    brochures :: Vector         # brochures from kp
    π :: Vector{Float64}        # hist productivity
    f :: Vector{Float64}        # hist market share
    μ :: Float64                # hist markup
    Π :: Vector{Float64}        # hist profits
    cI :: Float64               # internal funds for investments
    ΔDeb :: Float64             # changes in debt level
    balance :: Balance          # balance sheet
end


function initialize_cp(
    id :: Int, 
    machine_struct, 
    n_consrgood :: Int, 
    type_good :: String,
    n_init_emp_cp :: Int
    )

    cp = ConsumerGoodProducer(
        id,                     # global id
        type_good,              # type of good produced by producer
        [],                     # p: hist prices
        [],                     # c: hist cost
        [],                     # RD: hist R&D spending
        [1100],                 # D: hist demand
        1100,                   # Dᵉ exp demand
        Vector(),               # hh_queue: vector containing ids of demanding households           
        # [36e3],               # N: hist inventory
        40,                     # Nᵈ: desired inventory 
        [40],                   # Q: hist production
        40,                     # Qᵉ: exp production
        [],                     # I: hist investments
        [machine_struct],       # Ξ: capital stock
        40,                     # K: total amount of machines
        [],                     # RS: list of to-be replaced machines
        [],                     # employees: employees
        0,                      # L: labor units in company
        n_init_emp_cp * 100,    # Lᵉ: exp labor force
        0,                      # ΔLᵈ: desired change in labor force
        # 0,                      # P_FE: probability of getting fired while employed
        [1.0],                  # w: wage level
        1.0,                    # wᴼ: offered wage
        [],                     # brochures from kp
        [],                     # π: hist productivity
        [1/n_consrgood],        # f: market share
        0.05,                   # μ: hist markup
        [],                     # Π: hist profits
        0,                      # cI: internal funds for investments
        0,                      # ΔDeb: changes in debt level
        Balance(               
                1000,           # - N: inventory
                0.0,            # - K: capital
                0.0,            # - NW: liquid assets
                0.0,            # - Deb: debt
                0.0             # - EQ: equity
            )
    )
    return cp
end


"""
Plans production amounts for consumer good producer (short term)
    - updates ST expected demand
    - determines ST production goals
    - based on ST, set labor demand
"""
function plan_production_cp!(
    cp::ConsumerGoodProducer, 
    global_param,
    model::ABM
    )

    # update amount of owned capital
    update_K!(cp)

    # update expected demand
    cp.Dᵉ = global_param.ωD * cp.D[end] + (1 - global_param.ωD) * cp.Dᵉ

    # println(cp.Dᵉ)

    # determine desired short-term production
    Qˢ = cp.Dᵉ + cp.Nᵈ - cp.balance.N

    # println(cp.balance.N)

    # print(Qˢ)

    # update average productivity
    update_π!(cp)

    # compute corresponding change in labor stock
    # total_prod = sum(map(machine -> machine.A * machine.freq, cp.Ξ))
    total_K = sum(map(machine -> machine.A * machine.freq, cp.Ξ))
    # println(total_prod)
    # println("L:", cp.L)
    # println(Qˢ)
    # println(Qˢ/total_prod)

    # println(cp.balance.N, " ", cp.Dᵉ, " ", Qˢ/cp.π[end])

    # ΔLᵈ = Qˢ/cp.π[end] - cp.L
    ΔLᵈ = Qˢ/total_K - cp.L
    if ΔLᵈ < -cp.L
        cp.ΔLᵈ = -cp.L
    else
        cp.ΔLᵈ = ΔLᵈ
    end

    # Update average wage w̄
    update_wage_level_p!(cp, model)

    # Update markup μ
    compute_μ_cp!(cp, global_param.υ, global_param.μ1)

    # update cost of production c
    c = compute_c_cp!(cp, Qˢ)
    push!(cp.c, c)

    # Compute price
    p = (1 + cp.μ) * c
    push!(cp.p, p)
end


"""
Plans production amounts for consumer good producer (long term)
    - updates LT expected demand
    - updates LT labor supply 
    - determines LT production goals
    - based on LT, set investment amount
"""
function plan_investment_cp!(
    cp::ConsumerGoodProducer, 
    global_param, 
    all_kp::Vector{Int},
    model::ABM
    )

    # choose kp
    p_star, c_star, chosen_kp_id, Aᵈ = choose_producer_cp(cp, global_param.b, all_kp, model)

    # plan replacement investments
    plan_replacement_cp!(cp, global_param, p_star, c_star)

    # update LT demand
    cp.Qᵉ = global_param.ωQ * cp.Q[end] + (1 - global_param.ωQ) * cp.Qᵉ

    # println("Q ", cp.Qᵉ, " ", cp.Q[end])

    # update expected labor supply
    cp.Lᵉ = global_param.ωL * cp.L[end] + (1 - global_param.ωL) * cp.Lᵉ

    # println("L ", cp.Lᵉ, " ", cp.L[end])

    # compute desired capital stock expansion
    # Kᵈ = (cp.Qᵉ / cp.Lᵉ - sum(map(m -> m.A * m.freq, cp.Ξ)))/Aᵈ
    Kᵈ = (cp.Qᵉ / cp.Lᵉ - cp.π[end])/Aᵈ

    EIᵈ = Kᵈ - sum(map(machine -> machine.freq, cp.Ξ))
    
    Iₜ = EIᵈ + sum(map(machine -> machine.freq, cp.RS))

    # println(EIᵈ, " ", Iₜ, " ", EIᵈ)

    # TODO does not check if funds are available
    if Iₜ > 0
        order_machines_cp!(cp, chosen_kp_id, Iₜ, model)
    end

end


function update_K!(
    cp::ConsumerGoodProducer
    )
    cp.K = sum(map(machine -> machine.freq, cp.Ξ))
end


function update_π!(
    cp::ConsumerGoodProducer
    )
    π = sum(map(machine -> machine.freq / cp.K, cp.Ξ))
    push!(cp.π, π)
end


# Dosi et al (2013) Eq. 17, computes cost of production
cop(p_t, c_t, b) = p_t + b * c_t

"""
Plans replacement investment based on age machines and available new machines
"""
function plan_replacement_cp!(
    cp::ConsumerGoodProducer, 
    global_param,
    p_star::Float64, 
    c_star::Float64
    )

    cp.RS = []

    for machine in cp.Ξ
        if p_star/(machine.A - c_star) < global_param.b || machine.age >= global_param.η
            push!(cp.RS, machine)
        end
    end

end


function choose_producer_cp(
    cp::ConsumerGoodProducer, 
    b::Int, 
    all_kp::Vector{Int},
    model::ABM
    )

    if (length(cp.brochures) == 0)
        # in case of no brochures, pick a random kp
        chosen_kp_id = sample(all_kp)
        brochure = model[chosen_kp_id].brochure
        p_star = brochure[2]
        c_star = brochure[3]
        A_star = brochure[4]
        
        return p_star, c_star, chosen_kp_id, A_star
    end

    # take first producer as best, then compare to other producers
    chosen_kp_id = cp.brochures[1][1]
    p_star = cp.brochures[1][2]
    c_star = cp.brochures[1][3]
    A_star = cp.brochures[1][4]
    cop_star = cop(p_star, c_star, b)

    for brochure in cp.brochures

        p_h = brochure[2]
        c_h = brochure[3]
        A_h = brochure[4]
        potential_cop = cop(p_h, c_h, b)

        # if cheaper, choose cheaper option
        if potential_cop < cop_star
            chosen_kp_id = brochure[1]
            p_star = p_h
            c_star = c_h
            A_star = A_h
            cop_star = potential_cop
        end

    end

    return p_star, c_star, chosen_kp_id, A_star
end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(cp :: AbstractAgent)

    # compute weighted sum of machines
    # Ā = sum(map(machine -> machine.freq * machine.A, cp.Ξ))

    # compute total production amount
    Q = cp.π[end] * cp.L
    push!(cp.Q, Q)
    
    # change inventory, will be amount households can buy from
    cp.balance.N += Q
    # push!(cp.balance.N, N)

end


function compute_μ_cp!(
    cp::ConsumerGoodProducer, 
    υ::Float64, 
    μ1::Float64
    )

    if (length(cp.f) > 2)
        cp.μ = cp.μ * (1 + υ * (cp.f[end] - cp.f[end-1])/cp.f[end-1])
    else
        cp.μ = μ1
    end
end


function compute_c_cp!(
    cp::ConsumerGoodProducer, 
    Qˢ::Float64
    )::Float64

    total_K = sum(map(machine -> machine.freq, cp.Ξ))
    Ā = sum(map(m -> m.A * m.freq / total_K, cp.Ξ))
    # c = (cp.w̄[end] * Qˢ / Ā) / Qˢ
    c = cp.w̄[end] / cp.π[end]
    return c
end


function order_machines_cp!(
    cp::ConsumerGoodProducer,
    chosen_kp_id::Int, 
    Iₜ :: Float64,
    model::ABM
    )

    order = (cp.id, Iₜ)
    push!(model[chosen_kp_id].orders, order)

end


"""
Decides if enough inventory to send orders hh, 
    if yes: transact, if no: send back order unfulfilled or partially fulfilled
"""
function send_orders_cp!(
    cp::ConsumerGoodProducer,
    model::ABM
    )

    # Loop over orders in queue
    for order in cp.hh_queue

        hh_id = order[1]
        q = order[2]

        # Check if enough inventory available
        share_fulfilled = 0.0
        if cp.balance.N > q
            share_fulfilled = 1.0
            cp.D[end] += q
            cp.balance.N -= q
        elseif cp.balance.N > 0.0 && cp.balance.N < q
            share_fulfilled = cp.balance.N / q
            cp.D[end] += cp.balance.N
            cp.balance.N = 0.0
        end

        tot_price = (q * share_fulfilled) * cp.p[end]

        receive_order_hh!(model[hh_id], cp.id, tot_price, share_fulfilled)
    end
end


# """
# Determines if amount of demanded goods are available and transacts
# """
# function transact_cp!(
#     cp::ConsumerGoodProducer, 
#     hh::Household, 
#     N_G::Float64,
#     τˢ::Float64
#     )::Tuple{Bool, Float64}

#     # Check if enough inventory available, then transact
#     if N_G < cp.balance.N[end]

#         cp.balance.N -= N_G
#         cp.D[end] += N_G
#         hh.C -= cp.p[end] * (1 + τˢ) * N_G

#         sales_tax = cp.p[end] * τˢ * N_G

#         return true, sales_tax
#     else 
#         cp.D[end] += cp.balance.N
#         hh.C -= cp.p[end] * (1 + τˢ) * cp.balance.N
#         cp.balance.N = 0

#         sales_tax = cp.p[end] * τˢ * cp.balance.N
        
#         return false, sales_tax
#     end

# end


function compute_profits_cp!(
    cp::ConsumerGoodProducer
    )
    # TODO incorporate interest rates
    # println(cp.D[end])
    Π = cp.D[end] * cp.p[end] - cp.c[end] * cp.D[end]
    # println("Π ", Π)
    push!(cp.Π, Π)
end


function reset_brochures_cp!(
    cp::ConsumerGoodProducer
    )
    cp.brochures = []
end


function reset_queue_cp!(
    cp::ConsumerGoodProducer
    )
    cp.hh_queue = Vector()
end


"""
Lets cp receive machine and include it in capital stock.
"""
function receive_machines!(
    cp::ConsumerGoodProducer, 
    machine::Machine, 
    Iₜ::Float64,
    )

    # remove machines that were written off

    # println("yeet2 ", cp.Ξ)
    # println("yoot ", cp.RS)
    # println("yaat ", filter(m -> m ∉ cp.RS, cp.Ξ))

    # println("before ", length(cp.Ξ), " ", length(cp.RS))
    filter!(machine -> machine ∉ cp.RS, cp.Ξ)
    # println("after ", length(cp.Ξ))
    # println()

    # add new machine to capital stock
    push!(cp.Ξ, machine)

    # println("yeet3 ", cp.Ξ)

    cp.balance.NW -= Iₜ

    # TODO permutate the balance (in terms of capital and debt)

end