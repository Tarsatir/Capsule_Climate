"""
Defines struct for consumer good producer
"""
mutable struct ConsumerGoodProducer <: AbstractAgent
    id :: Int                       # id
    type_good :: String             # type of good produced by producer 
    p :: Vector{Float64}            # hist prices
    c :: Vector{Float64}            # hist cost
    RD :: Vector{Float64}           # R&D spending
    D :: Vector{Float64}            # hist demand
    Dᵉ :: Float64                   # exp demand
    hh_queue :: Vector              # vector containing orders of demanding households
    Nᵈ :: Float64                   # desired inventory
    N_goods :: Float64              # inventory in good units
    Q :: Vector{Float64}            # hist production
    Qᵉ :: Float64                   # exp production
    Qˢ :: Float64                   # desired short-term production

    # Investments
    Iᵈ :: Float64                   # desired total investments
    EIᵈ :: Float64                  # desired expansionary investments
    RSᵈ :: Float64                  # desired replacement investments
    RS :: Vector{Machine}           # list of to-be replaced machines
    chosen_kp_id :: Int             # id of chosen kp
    Deb_installments :: Vector{Float64}   # installments of debt repayments

    Ξ :: Vector{Machine}            # capital stock
    employees:: Vector{Int}         # employees list
    L :: Float64                    # labor units
    Lᵉ:: Float64                    # exp labor force
    ΔLᵈ :: Float64                  # desired change in labor force
    w̄ :: Vector{Float64}            # wage level
    wᴼ :: Float64                   # offered wage
    wᴼₑ :: Float64                  # expected offered wage 
    brochures :: Vector             # brochures from kp
    π :: Vector{Float64}            # hist productivity
    f :: Vector{Float64}            # hist market share
    μ :: Float64                    # markup rate
    Π :: Vector{Float64}            # hist profits
    cI :: Float64                   # internal funds for investments
    balance :: Balance              # balance sheet
    curracc :: FirmCurrentAccount   # current account
end


function initialize_cp(
    id :: Int, 
    machines::Vector{Machine}, 
    n_consrgood :: Int, 
    type_good :: String,
    n_init_emp_cp :: Int
    )

    balance = Balance(               
        0,                          # - N: inventory (in money units)
        n_init_emp_cp * 100,        # - K: capital
        1000.0,                     # - NW: liquid assets
        0.0,                        # - Deb: debt
        0.0                         # - EQ: equity
    )

    curracc = FirmCurrentAccount(0,0,0,0,0,0)

    cp = ConsumerGoodProducer(
        id,                         # global id
        type_good,                  # type of good produced by producer
        [],                         # p: hist prices
        [],                         # c: hist cost
        [],                         # RD: hist R&D spending
        [1100],                     # D: hist demand
        1100,                       # Dᵉ exp demand
        Vector(),                   # hh_queue: vector containing ids of demanding households           
        1000,                       # Nᵈ: desired inventory
        1000,                       # N_goods: inventory in goods 
        [1100],                     # Q: hist production
        1100,                       # Qᵉ: exp production
        0,                          # Qˢ: desired short-term production

        0,                          # Iᵈ: desired investments
        0,                          # EIᵈ: desired expansionary investments
        0,                          # RSᵈ: desired replacement investments
        [],                         # RS: list of to-be replaced machines
        0,                          # chosen_kp_id: id of chosen kp
        Vector{Float64}([0.0,       # Deb_installments: vector containing debt repayments
                         0.0, 
                         0.0]),

        machines,                   # Ξ: capital stock
        [],                         # employees: employees
        0,                          # L: labor units in company
        n_init_emp_cp * 100,        # Lᵉ: exp labor force
        0,                          # ΔLᵈ: desired change in labor force
        [1.0],                      # w: wage level
        1.0,                        # wᴼ: offered wage
        1.0,                        # wᴼₑ: expected offered wage
        [],                         # brochures from kp
        [],                         # π: hist productivity
        [0.5/n_consrgood],          # f: market share
        0.05,                       # μ: hist markup
        [],                         # Π: hist profits
        0,                          # cI: internal funds for investments
        balance,                    # balance
        curracc                     # current account
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

    # Update amount of owned capital
    update_K!(cp)

    # Update expected demand
    cp.Dᵉ = global_param.ωD * cp.D[end] + (1 - global_param.ωD) * cp.Dᵉ

    # println(cp.Dᵉ)

    # Determine desired short-term production
    cp.Qˢ = cp.Dᵉ + cp.Nᵈ - cp.N_goods

    # Update average productivity
    update_π!(cp)

    # Compute corresponding change in labor stock
    compute_ΔLᵈ_cp!(cp)

    # Update average wage w̄
    update_wage_level_p!(cp, model)

    # Update markup μ
    compute_μ_cp!(cp, global_param.υ, global_param.μ1)

    # Update cost of production c
    compute_c_cp!(cp)

    # Compute price
    compute_p_cp!(cp)
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
    global_param::GlobalParam, 
    all_kp::Vector{Int},
    model::ABM
    )

    # Choose kp
    p_star, Aᵈ = choose_producer_cp(cp, global_param.b, all_kp, model)

    # Plan replacement investments
    plan_replacement_cp!(cp, global_param, p_star, Aᵈ)

    # Update LT demand
    cp.Qᵉ = global_param.ωQ * cp.Q[end] + (1 - global_param.ωQ) * cp.Qᵉ

    # Update expected labor supply
    cp.Lᵉ = global_param.ωL * cp.L[end] + (1 - global_param.ωL) * cp.Lᵉ

    # Compute desired capital stock expansion
    Kᵈ = cp.Qᵉ / cp.π[end]

    cp.EIᵈ = max(Kᵈ - cp.balance.K, 0)
    cp.RSᵈ = max(sum(map(machine -> machine.freq, cp.RS)), 0)
    cp.Iᵈ = max(cp.EIᵈ + cp.RSᵈ, 0)
end


"""
Makes funding allocation based on desired production and investment. 
If not enough funding available, changes desired production and investment
    to a possible level.
"""
function funding_allocation_cp!(
    cp::ConsumerGoodProducer,
    Λ::Float64,
    model::ABM
    )

    # Compute how much debt cp can make additionally
    poss_add_Deb = max(cp.D[end] * cp.p[end] * Λ - cp.balance.Deb, 0)

    # println(cp.p[end])

    # Compute the expected total cost of desired production
    # TODO: move this to hiring desicion (change n labor in hh)
    ΔL_ceil = 100 * ceil(cp.ΔLᵈ / 100)
    total_cop = max(cp.L * cp.w̄[end] + ΔL_ceil * cp.wᴼₑ, 0)

    # println(cp.balance.NW, " ", total_cop, " ", cp.Iᵈ, " ", poss_add_Deb)
    # println(cp.balance.Deb)

    # Check if enough internal funds available
    if cp.balance.NW > total_cop && cp.Iᵈ == 0
        # Production can be funded from NW, no investment desired
        # Already available in cash amount, no need to borrow

        # Use remaining NW to pay back debts
        # payback = min(cp.balance.Deb, cp.balance.NW - total_cop)
        # payback_debt_p!(cp, payback)

    elseif cp.balance.NW > total_cop + cp.Iᵈ
        # Production and full investment can be financed from NW
        # Already available in cash amount, no need to borrow
        
        # Use remaining NW to pay back debts
        # payback = min(cp.balance.Deb, cp.balance.NW - total_cop)
        # payback_debt_p!(cp, payback)

    elseif cp.balance.NW + poss_add_Deb > total_cop + cp.Iᵈ
        # Production and full investment partially funded from Deb
        # Compute required extra cash and borrow
        req_cash = total_cop + cp.Iᵈ - cp.balance.NW
        # borrow_funds_p!(cp, req_cash)
    else
        # Production or investment constrained by financial means
        if cp.balance.NW + poss_add_Deb < total_cop
            # Production constrained. Decrease production target and labor demand,
            # No investments possible

            # Recompute labor that can be hired
            cp.ΔLᵈ = (cp.balance.NW + poss_add_Deb - cp.w̄[end] * cp.L)/cp.wᴼₑ

            # borrow_funds_p!(cp, poss_add_Deb)

            # Set all desired investments to zero
            cp.Iᵈ = 0
            cp.EIᵈ = 0
            cp.RSᵈ = 0
            cp.RS = Vector{Machine}()

        else
            # Production possible from NW and Deb, investment constrained

            # Compute extra required cash for production, borrow amount
            req_cash = total_cop - cp.balance.NW

            if cp.balance.NW + poss_add_Deb < total_cop + cp.EIᵈ
                # Partial expansion investment, no replacement
                poss_EI = cp.balance.NW + poss_add_Deb - total_cop - req_cash
                req_cash += poss_EI

                # Update desired investments
                cp.Iᵈ = poss_EI
                cp.EIᵈ = poss_EI
                cp.RSᵈ = 0
                cp.RS = Vector{Machine}()
            else
                # Full expansion investment, partial replacement
                # poss_EI = cp.balance.NW + poss_add_Deb - total_cop - req_cash
                # poss_RS = cp.balance.NW + poss_add_Deb - total_cop - req_cash - poss_EI

                # P_RS = model[cp.chosen_kp_id].p[end]
                # TODO let cp replace partial capital goods

                req_cash += cp.EIᵈ

                cp.Iᵈ = cp.EIᵈ
                cp.RSᵈ = 0
                cp.RS = Vector{Machine}()

            end

            # borrow_funds_p!(cp, req_cash)
        end
    end
end


function update_K!(
    cp::ConsumerGoodProducer
    )

    cp.balance.K = sum(map(machine -> machine.freq, cp.Ξ))
end


function update_π!(
    cp::ConsumerGoodProducer
    )

    π = sum(map(machine -> (machine.freq * machine.A) / cp.balance.K, cp.Ξ))
    push!(cp.π, π)
end


"""
Adaptively updates expected offered wage
"""
function update_wᴼₑ(
    cp::ConsumerGoodProducer,
    ωW::Float64
    )

    cp.wᴼₑ = ωW * cp.wᴼₑ + (1 - ωW) * cp.wᴼ
end


# Dosi et al (2013) Eq. 17, computes cost of production
cop(p_t, c_t, b) = p_t + b * c_t


"""
Plans replacement investment based on age machines and available new machines
"""
function plan_replacement_cp!(
    cp::ConsumerGoodProducer, 
    global_param::GlobalParam,
    p_star::Float64, 
    Aᵈ::Float64
    )

    cp.RS = Vector{Machine}()

    for machine in cp.Ξ
        c_A = cp.w̄[end] / machine.A
        c_star = cp.w̄[end] / Aᵈ
        if c_A > c_star
            if p_star/(c_A - c_star) <= global_param.b || machine.age >= global_param.η
                push!(cp.RS, machine)
            end
        end
    end
end


function choose_producer_cp(
    cp::ConsumerGoodProducer, 
    b::Int, 
    all_kp::Vector{Int},
    model::ABM
    )::Tuple{Float64, Float64}

    if (length(cp.brochures) == 0)
        # in case of no brochures, pick a random kp
        cp.chosen_kp_id = sample(all_kp)
        brochure = model[cp.chosen_kp_id].brochure
        p_star = brochure[2]
        c_star = brochure[3]
        A_star = brochure[4]
        
        return p_star, A_star
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

    cp.chosen_kp_id = chosen_kp_id

    return p_star, A_star
end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(
    cp::ConsumerGoodProducer
    )

    # Compute total production amount
    Q = min(cp.π[end] * cp.L, cp.balance.K, cp.Qˢ)
    push!(cp.Q, Q)
    
    # Change inventory, will be amount households can buy from
    cp.N_goods += Q
end


function compute_μ_cp!(
    cp::ConsumerGoodProducer, 
    υ::Float64, 
    μ1::Float64
    )

    # First check if market share is nonzero, else liquidate firm
    if cp.f[end] <= 0.001
        # TODO: liquidate firm
        cp.f[end] = 0.001
    end

    # Compute market share
    if length(cp.f) > 2
        cp.μ = cp.μ * (1 + υ * (cp.f[end] - cp.f[end-1])/cp.f[end-1])
    else
        cp.μ = μ1
    end
end


function compute_c_cp!(
    cp::ConsumerGoodProducer,
    )

    c = cp.w̄[end] / cp.π[end]
    push!(cp.c, c)
end


function compute_p_cp!(
    cp::ConsumerGoodProducer
    )

    p = (1 + cp.μ) * cp.c[end]
    push!(cp.p, p)
end


function compute_ΔLᵈ_cp!(
    cp::ConsumerGoodProducer
    )

    L = 100 # TODO: make this a global parameter

    ΔLᵈ = cp.Qˢ/cp.π[end] - cp.L

    if ΔLᵈ < -cp.L
        cp.ΔLᵈ = -cp.L
    else
        cp.ΔLᵈ = L * ceil(ΔLᵈ / L)
    end
end


function order_machines_cp!(
    cp::ConsumerGoodProducer,
    model::ABM
    )

    if cp.Iᵈ > 0
        order = (cp.id, cp.Iᵈ)
        push!(model[cp.chosen_kp_id].orders, order)
    end
end


"""
Decides if enough inventory to send orders hh, 
    if yes: transact, if no: send back order unfulfilled or partially fulfilled
"""
function send_orders_cp!(
    cp::ConsumerGoodProducer,
    model::ABM
    )

    # all_q = sum(map(o -> o[2], cp.hh_queue))
    # avg_C = mean(map(o -> model[o[1]].C[end], cp.hh_queue))
    # println(cp.Q[end], " ", all_q, " ", avg_C)

    push!(cp.D, 0)
    cp.curracc.S = 0

    # Loop over orders in queue
    for order in cp.hh_queue

        hh_id = order[1]
        q = order[2]

        # Check if enough inventory available
        share_fulfilled = 0.0
        if cp.N_goods > q
            share_fulfilled = 1.0
            cp.D[end] += q
            cp.N_goods -= q
        elseif cp.N_goods > 0.0 && cp.N_goods < q
            share_fulfilled = cp.N_goods / q
            cp.D[end] += cp.N_goods
            cp.N_goods = 0.0
        end

        tot_price = (q * share_fulfilled) * cp.p[end]
        cp.curracc.S += tot_price
        receive_order_hh!(model[hh_id], cp.id, tot_price, share_fulfilled)

        if cp.N_goods == 0.0
            return
        end
    end
end


function compute_profits_cp!(
    cp::ConsumerGoodProducer,
    r::Float64
    )

    # Compute total cost of investments
    cp.curracc.int_Deb = cp.balance.Deb * r
 
    Π = (cp.p[end] - cp.c[end]) * cp.D[end] - cp.curracc.int_Deb
    push!(cp.Π, Π)
end


"""
Updates desired inventory to be a share of the previous demand
"""
function update_Nᵈ_cp!(
    cp::ConsumerGoodProducer,
    Nᵈ_share::Float64
    )

    cp.Nᵈ = Nᵈ_share * cp.Dᵉ
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


function increase_machine_age_cp!(
    cp::ConsumerGoodProducer
    )

    for machine in cp.Ξ
        machine.age += 1
    end
end


"""
Lets cp receive machine and include it in capital stock.
"""
function receive_machines!(
    cp::ConsumerGoodProducer, 
    machine::Machine, 
    Iₜ::Float64,
    )

    # Remove machines that were written off
    # println("1 ", length(cp.Ξ), " ", length(cp.RS))
    filter!(machine -> machine ∉ cp.RS, cp.Ξ)
    # println("2 ", length(cp.Ξ))

    # Add new machine to capital stock
    push!(cp.Ξ, machine)
    # cp.balance.NW -= Iₜ

    # Add debt to future installments
    cp.Deb_installments[1] += Iₜ / 3
    cp.Deb_installments[2] += Iₜ / 3
    cp.Deb_installments[3] += Iₜ / 3


    # TODO permutate the balance (in terms of capital and debt)

end


"""
Repays outstanding debts at start of period, shifts debt installments
"""
function repay_debt_installments_cp!(
    cp::ConsumerGoodProducer
    )

    to_be_paid = cp.Deb_installments[1]
    cp.curracc.rep_Deb += to_be_paid

    cp.Deb_installments[1] = cp.Deb_installments[2]
    cp.Deb_installments[2] = cp.Deb_installments[3]
    cp.Deb_installments[3] = 0.0
end