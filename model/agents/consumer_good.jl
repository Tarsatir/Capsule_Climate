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
    order_queue :: Vector           # vector containing orders of demanding households
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

    Ξ :: Vector{Machine}            # machines
    n_machines :: Float64                    # total freq of machines # TODO rename
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
    n_init_emp_cp :: Int,
    μ1::Float64,
    ι::Float64
    )

    balance = Balance(               
        0,                          # - N: inventory (in money units)
        n_init_emp_cp * 100,        # - K: capital
        100.0,                      # - NW: liquid assets
        0.0,                        # - Deb: debt
        0.0                         # - EQ: equity
    )

    curracc = FirmCurrentAccount(0,0,0,0,0,0)

    cp = ConsumerGoodProducer(
        id,                         # global id
        type_good,                  # type of good produced by producer
        [1.0 * (1 + μ1)],           # p: hist prices
        [1.0],                      # c: hist cost
        [],                         # RD: hist R&D spending
        [1100],                     # D: hist demand
        1100,                       # Dᵉ exp demand
        Vector(),                   # hh_queue: vector containing ids of demanding households           
        ι * 1100,                   # Nᵈ: desired inventory
        0,                          # N_goods: inventory in goods 
        [1210],               # Q: hist production
        1210,                 # Qᵉ: exp production
        0,                          # Qˢ: desired short-term production

        0,                          # Iᵈ: desired investments
        0,                          # EIᵈ: desired expansionary investments
        0,                          # RSᵈ: desired replacement investments
        [],                         # RS: list of to-be replaced machines
        0,                          # chosen_kp_id: id of chosen kp
        Vector{Float64}([0.0,       # Deb_installments: vector containing debt repayments
                         0.0, 
                         0.0]),

        machines,                   # Ξ: machines
        0,                          # n_machines: total number of machines
        [],                         # employees: employees
        n_init_emp_cp * 100,        # L: labor units in company
        n_init_emp_cp * 100,        # Lᵉ: exp labor force
        0,                          # ΔLᵈ: desired change in labor force
        [1.0],                      # w: wage level
        1.0,                        # wᴼ: offered wage
        1.0,                        # wᴼₑ: expected offered wage
        [],                         # brochures from kp
        [],                         # π: hist productivity
        [2/n_consrgood],            # f: market share
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
    global_param::GlobalParam,
    model::ABM
    )

    # Update amount of owned capital and desired inventories
    update_n_machines_cp!(cp)
    update_Nᵈ_cp!(cp, global_param.ι)

    # Compute expected demand
    update_Dᵉ_cp!(cp, global_param.ωD)

    # Compute desired short-term production
    update_Qˢ_cp!(cp)
    # println(cp.Qˢ)

    # Update average productivity
    compute_π_cp!(cp)

    # Compute corresponding change in labor stock
    update_ΔLᵈ_cp!(cp)

    # Update average wage w̄
    update_w̄_p!(cp, model)

    # Update markup μ
    update_μ_cp!(cp, global_param.υ, global_param.μ1)

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
    # p_star, Aᵈ = choose_producer_cp!(cp, global_param.b, all_kp, model)
    brochure = choose_producer_cp!(cp, global_param.b, all_kp, model)

    # Plan replacement investments
    plan_replacement_cp!(cp, global_param, brochure)

    # Update LT production
    update_Qᵉ_cp!(cp)

    # Update expected labor supply
    # update_Lᵉ_cp!(cp)
    # TODO: see if I still need this for something

    # Compute desired capital stock expansion
    Kᵈ = cp.Qᵉ / cp.π[end]
    cp.EIᵈ = min(max(Kᵈ - cp.n_machines, 0), cp.n_machines * global_param.Kg_max)
    cp.RSᵈ = max(sum(map(machine -> machine.freq, cp.RS)), 0)
    # println(Kᵈ, " ", cp.EIᵈ, " ", cp.RSᵈ, " ", cp.n_machines)
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
    ωW::Float64,
    model::ABM
    )

    # Update offered wage expectation
    update_wᴼₑ_cp!(cp, ωW)

    # Compute how much debt cp can make additionally
    poss_add_Deb = max((cp.D[end] * cp.p[end - 1]) * Λ - sum(cp.Deb_installments), 0)

    # Compute the expected total cost of desired production
    total_cop = cp.L * cp.w̄[end] + cp.ΔLᵈ * cp.wᴼₑ

    # Check if enough internal funds available
    if cp.balance.NW > total_cop + cp.Iᵈ
        # Production and full investment can be financed from NW
        # Already available in cash amount, no need to borrow
        cp.cI = cp.Iᵈ

    elseif cp.balance.NW + poss_add_Deb > total_cop + cp.Iᵈ
        # Production and full investment partially funded from Deb
        # Compute required extra cash and borrow
        # req_cash = total_cop + cp.Iᵈ - cp.balance.NW
        cp.cI = max(cp.balance.NW - total_cop, 0)

    else
        # Production or investment constrained by financial means
        if cp.balance.NW + poss_add_Deb < total_cop
            # Production constrained. Decrease production target and labor demand,
            # No investments possible

            # Recompute labor that can be hired
            if cp.balance.NW + poss_add_Deb > cp.L * cp.w̄[end]
                # Current stock of labor has to be brought down, no additional hires
                cp.ΔLᵈ = (cp.balance.NW + poss_add_Deb) / cp.w̄[end] - cp.L
            else
                # Desired increase in labor stock has to be brought down
                cp.ΔLᵈ = (cp.balance.NW + poss_add_Deb - cp.w̄[end] * cp.L)/cp.wᴼₑ
            end

            # Set all desired investments to zero
            cp.Iᵈ = 0
            cp.EIᵈ = 0
            cp.RSᵈ = 0
            cp.RS = Vector{Machine}()

            cp.cI = 0.0

        else
            # Production possible from NW and Deb, investment constrained

            # Compute extra required cash for production, borrow amount
            # req_cash = total_cop - cp.balance.NW

            if cp.balance.NW + poss_add_Deb < total_cop + cp.EIᵈ
                # Partial expansion investment, no replacement
                poss_EI = cp.balance.NW + poss_add_Deb - total_cop
                # req_cash += poss_EI

                # Update desired investments
                cp.Iᵈ = poss_EI
                cp.EIᵈ = poss_EI
                cp.RSᵈ = 0
                cp.RS = Vector{Machine}()

                cp.cI = max(cp.balance.NW - total_cop, 0)
            else
                # Full expansion investment, partial replacement
                # poss_EI = cp.balance.NW + poss_add_Deb - total_cop - req_cash
                # poss_RS = cp.balance.NW + poss_add_Deb - total_cop - req_cash - poss_EI

                # P_RS = model[cp.chosen_kp_id].p[end]
                # TODO let cp replace partial capital goods

                # req_cash += cp.EIᵈ

                cp.Iᵈ = cp.EIᵈ
                cp.RSᵈ = 0
                cp.RS = Vector{Machine}()

                cp.cI = max(cp.balance.NW - total_cop, 0)

            end
        end
    end
end


"""
Updates expected demand based Dᵉ
"""
function update_Dᵉ_cp!(
    cp::ConsumerGoodProducer,
    ωD::Float64
    )

    if length(cp.D) > 2
        cp.Dᵉ = ωD * cp.Dᵉ + (1 - ωD) *  (2*cp.D[end] - cp.D[end-1])
    else
        cp.Dᵉ = ωD * cp.Dᵉ + (1 - ωD) *  cp.D[end]
    end
end


"""
Updates desired short-term production Qˢ
"""
function update_Qˢ_cp!(
    cp::ConsumerGoodProducer
    )

    cp.Qˢ = cp.Dᵉ + cp.Nᵈ - cp.N_goods
    # if length(cp.D) > 1
        # println(cp.Qˢ, " ", cp.Dᵉ, " ", cp.Nᵈ, " ", cp.N_goods, " ", cp.D[end], " ", cp.D[end-1])
    # end
    # println(cp.Dᵉ, " ", cp.N_goods)
end


"""
Updates expected long-term production Qᵉ
"""
function update_Qᵉ_cp!(
    cp::ConsumerGoodProducer
    )

    # TODO: very large sensitivity to this parameter: check out!
    
    # if length(cp.D) > 2
    #     Qg = cp.Q[end] * (1 + (cp.D[end] - cp.D[end-1]) / cp.D[end-1])
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * Qg
    # else
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end]
    # end

    if length(cp.Π) > 2 && cp.Π[end-1] != 0
        Qg = cp.Q[end] * (1 + (cp.Π[end] - cp.Π[end-1]) / cp.Π[end-1])
        cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * Qg
    else
        cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end]
    end
end


# """
# Updates expected long-term labor supply Lᵉ
# """
# function update_Lᵉ_cp!(
#     cp::ConsumerGoodProducer
#     )

#     cp.Lᵉ = global_param.ωL * cp.L[end] + (1 - global_param.ωL) * cp.Lᵉ
# end


"""
Updates capital stock n_machines
"""
function update_n_machines_cp!(
    cp::ConsumerGoodProducer
    )

    cp.n_machines = sum(map(machine -> machine.freq, cp.Ξ))
end


"""
Updates weighted producivity of machine stock π
"""
function compute_π_cp!(
    cp::ConsumerGoodProducer
    )

    # total_freq = sum(map(machine -> machine.freq, cp.Ξ))
    π = sum(map(machine -> (machine.freq * machine.A) / cp.n_machines, cp.Ξ))
    push!(cp.π, π)
end


"""
Adaptively updates expected offered wage
"""
function update_wᴼₑ_cp!(
    cp::ConsumerGoodProducer,
    ωW::Float64
    )

    cp.wᴼₑ = ωW * cp.wᴼₑ + (1 - ωW) * cp.wᴼ
end


"""
Plans replacement investment based on age machines and available new machines
"""
function plan_replacement_cp!(
    cp::ConsumerGoodProducer, 
    global_param::GlobalParam,
    # p_star::Float64, 
    # Aᵈ::Float64
    brochure
    )

    p_star = brochure[2]
    Aᵈ = brochure[4]

    cp.RS = Vector{Machine}()

    for machine in cp.Ξ
        c_A = cp.w̄[end] / machine.A
        c_star = cp.w̄[end] / Aᵈ
        if c_A > c_star
            if p_star/(c_A - c_star) <= global_param.b || machine.age >= global_param.η
            # if p_star/(c_A - c_star) <= global_param.b
                push!(cp.RS, machine)
            end
        end
    end
end


"""
Lets cp make decision for kp out of available kp in brochures.
"""
function choose_producer_cp!(
    cp::ConsumerGoodProducer, 
    b::Int, 
    all_kp::Vector{Int},
    model::ABM
    )

    # In case of no brochures, pick a random kp
    if (length(cp.brochures) == 0)
        cp.chosen_kp_id = sample(all_kp)
        brochure = model[cp.chosen_kp_id].brochure
        # p_star = brochure[2]
        # c_star = brochure[3]
        # A_star = brochure[4]
        
        # return p_star, A_star
        return brochure
    end

    # Take first producer as best, then compare to other producers
    chosen_kp_id = cp.brochures[1][1]
    p_star = cp.brochures[1][2]
    c_star = cp.brochures[1][3]
    A_star = cp.brochures[1][4]
    brochure_star = cp.brochures[1]
    cop_star = p_star + b * c_star

    for brochure in cp.brochures

        p_h = brochure[2]
        c_h = brochure[3]
        A_h = brochure[4]
        potential_cop = p_h + b * c_h

        # if cheaper, choose cheaper option
        if potential_cop < cop_star
            brochure_star = brochure
            chosen_kp_id = brochure[1]
            p_star = p_h
            c_star = c_h
            A_star = A_h
            cop_star = potential_cop
        end

    end

    cp.chosen_kp_id = chosen_kp_id

    # return p_star, A_star
    return brochure_star
end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(
    cp::ConsumerGoodProducer
    )

    # Compute total production amount
    Q = min(cp.π[end] * cp.L, cp.π[end] * cp.n_machines, cp.Qˢ)
    # println("L ", cp.π[end] * cp.L, " ", cp.n_machines, " ", cp.Qˢ)
    push!(cp.Q, Q)
    
    # Change inventory, will be amount households can buy from
    cp.N_goods += Q
end


"""
Computes the markup rate μ based on the market share f.
"""
function update_μ_cp!(
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


"""
Compute production cost per unit c
"""
function compute_c_cp!(
    cp::ConsumerGoodProducer,
    )

    # Check if any employees available
    if cp.L + cp.ΔLᵈ > 0
        # wᵉ = (cp.w̄[end] * cp.L + cp.wᴼₑ * cp.ΔLᵈ) / (cp.L + cp.ΔLᵈ)
        # println("c ", wᵉ)
        # c = wᵉ / cp.π[end]
        c = cp.w̄[end] / cp.π[end]
        push!(cp.c, c)
    end
end


"""
Computes price based on cost c and markup μ
"""
function compute_p_cp!(
    cp::ConsumerGoodProducer
    )

    p = (1 + cp.μ) * cp.c[end]
    # println("p ", p, " c ", cp.c[end], " μ ", cp.μ)
    push!(cp.p, p)
end


"""
Computes the desired change in labor supply ΔLᵈ
"""
function update_ΔLᵈ_cp!(
    cp::ConsumerGoodProducer
    )

    L = 100 # TODO: make this a global parameter

    ΔLᵈ = min(cp.Qˢ, cp.n_machines * cp.π[end])/cp.π[end] - cp.L

    if ΔLᵈ < -cp.L
        cp.ΔLᵈ = -cp.L
    else
        cp.ΔLᵈ = L * floor(ΔLᵈ / L)
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

    # Compute total demand for goods, regardless of if it can be fulfilled
    # Allows producers to produce more if production is insufficient
    total_D = sum(map(order -> order[2], cp.order_queue))
    push!(cp.D, total_D)

    # Loop over orders in queue, add realized sales S
    for order in cp.order_queue

        hh_id = order[1]
        q = order[2]

        if q > 0

            # Check if enough inventory available
            share_fulfilled = 0.0
            if cp.N_goods > q
                share_fulfilled = 1.0
                cp.N_goods -= q
            elseif cp.N_goods > 0.0 && cp.N_goods < q
                share_fulfilled = cp.N_goods / q
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
end


"""
Computes profit of cp in previous time period
"""
function compute_Π_cp!(
    cp::ConsumerGoodProducer,
    r::Float64,
    writeoffs=0::Float64
    )

    # Compute total cost of investments
    cp.curracc.int_Deb = cp.balance.Deb * r
    cp.curracc.TCL = cp.L * cp.w̄[end]

    # println(cp.curracc.S, " ", cp.curracc.TCL, " ", cp.curracc.int_Deb)
 
    Π = cp.curracc.S + cp.curracc.rev_dep - cp.curracc.TCL - cp.curracc.int_Deb
    # Π = cp.curracc.S + cp.curracc.rev_dep - cp.curracc.TCL - cp.curracc.int_Deb - writeoffs
    push!(cp.Π, Π)
end


"""
Updates desired inventory to be a share of the previous demand
"""
function update_Nᵈ_cp!(
    cp::ConsumerGoodProducer,
    ι::Float64
    )

    cp.Nᵈ = ι * cp.D[end]
end


function reset_brochures_cp!(
    cp::ConsumerGoodProducer
    )

    cp.brochures = []
end


function reset_queue_cp!(
    cp::ConsumerGoodProducer
    )

    cp.order_queue = Vector()
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
    # cp.curracc.TCI += Iₜ

    # Add debt to future installments
    add_debt = max(Iₜ - cp.cI, 0)
    borrow_funds_p!(cp, add_debt)

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