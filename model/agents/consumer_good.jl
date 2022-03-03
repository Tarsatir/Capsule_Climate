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
    mach_tb_repl :: Vector{Machine}           # list of to-be replaced machines
    chosen_kp_id :: Int             # id of chosen kp
    debt_installments :: Vector{Float64}   # installments of debt repayments

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
    π :: Float64                    # hist productivity
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
        0.0,                        # - debt: debt
        0.0                         # - EQ: equity
    )

    curracc = FirmCurrentAccount(0,0,0,0,0,0,0)

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
        [1210],                     # Q: hist production
        1210,                       # Qᵉ: exp production
        0,                          # Qˢ: desired short-term production
        0,                          # Iᵈ: desired investments
        0,                          # EIᵈ: desired expansionary investments
        0,                          # RSᵈ: desired replacement investments
        [],                         # mach_tb_repl: list of to-be replaced machines
        0,                          # chosen_kp_id: id of chosen kp
        fill(0.0, 4),               # debt installments
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
        1.0,                        # π: hist productivity
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

    # Update average productivity
    update_π_cp!(cp)

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
    brochure = choose_producer_cp!(cp, global_param.b, all_kp, model)

    # Plan replacement investments
    plan_replacement_cp!(cp, global_param, brochure)

    # Update LT production
    update_Qᵉ_cp!(cp)

    # Update expected labor supply
    # update_Lᵉ_cp!(cp)
    # TODO: see if I still need this for something

    # Compute desired capital stock expansion
    desired_n_machines = cp.Qᵉ / cp.π
    p_machine  = brochure[2]
    cp.EIᵈ = p_machine * min(max(desired_n_machines - cp.n_machines, 0), 
                             cp.n_machines * global_param.Kg_max)
    cp.RSᵈ = p_machine * max(sum(map(machine -> machine.freq, cp.mach_tb_repl)), 0)
    # println(Kᵈ, " ", cp.EIᵈ, " ", cp.RSᵈ, " ", cp.n_machines)
    cp.Iᵈ = max(cp.EIᵈ + cp.RSᵈ, 0)
end


"""
Checks funding restructions based on expected revenue and expenses. If not enough
    funding available in firm, decrease desired production, hiring or investments.
"""
function check_funding_restrictions_cp!(
    cp::ConsumerGoodProducer,
    Λ::Float64,
    r::Float64,
    ωW::Float64,
    model::ABM
    )

    # Update expected TCL
    TCLᵉ = (cp.L + cp.ΔLᵈ) * cp.w̄[end]
    if cp.ΔLᵈ > 0
        TCLᵉ = cp.w̄[end] * cp.L + cp.ΔLᵈ * cp.wᴼ
    end

    max_debt = Λ * cp.D[end] * cp.p[end - 1]
    max_add_debt = max_debt - cp.balance.debt

    # Check if cost of labor and investment can be financed from liquid assets
    NW_no_prod = cp.balance.NW + cp.Dᵉ * cp.p[end] + cp.curracc.rev_dep - cp.debt_installments[1] - cp.balance.debt * r
    
    if NW_no_prod > TCLᵉ + cp.Iᵈ
        # All cost of labor and investments can be paid from liquid assets. No additional debt needed.
        cp.cI = cp.Iᵈ

    elseif NW_no_prod > TCLᵉ && NW_no_prod - TCLᵉ < cp.Iᵈ
        # Cost of labor can be paid from liquid assets, investment has to be partially funded from debt.
        add_debt = cp.Iᵈ - (NW_no_prod - TCLᵉ)

        if add_debt > max_add_debt
            # If additional debt more than max additional debt, decrease investment amount
            cp.cI = 0.0

            if add_debt > cp.EIᵈ
                # Decrease amount of expansionary investment.
                cp.Iᵈ = max_add_debt
                cp.RSᵈ = 0.0
                cp.mach_tb_repl = Vector{Machine}()
            else
                # No replacement investment, only expansionary investment.
                cp.Iᵈ = cp.EIᵈ
                cp.RSᵈ = 0.0
                cp.mach_tb_repl = Vector{Machine}()
            end

        else
            # Additional debt less than maximum additional debt. Borrow to finance investment.
            cp.cI = cp.Iᵈ - add_debt
        end

    else
        # Cost of labor exceeds liquid assets. Check if enough additional debt available.
        # All investment cancelled.
        cp.cI = 0.0
        cp.Iᵈ = 0.0
        cp.EIᵈ = 0.0
        cp.RSᵈ = 0.0
        cp.mach_tb_repl = Vector{Machine}()

        if NW_no_prod + max_add_debt < TCLᵉ
            # Cost of labor exceeds expected liquid assets plus max additional debt. 
            # Decrease production quantity.
            poss_TCL = NW_no_prod + max_add_debt
            poss_L = poss_TCL / cp.w̄[end]
            cp.ΔLᵈ = poss_L - cp.L
            cp.Qˢ = min((cp.L + cp.ΔLᵈ) * cp.π[end], cp.n_machines * cp.π[end])
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
        cp.Dᵉ = ωD * cp.Dᵉ + (1 - ωD) *  (2 * cp.D[end] - cp.D[end-1])
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

    cp.Qˢ = max(cp.Dᵉ + cp.Nᵈ - cp.N_goods, 0)
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


"""
Updates expected long-term labor supply Lᵉ
"""
function update_Lᵉ_cp!(
    cp::ConsumerGoodProducer
    )

    cp.Lᵉ = global_param.ωL * cp.L[end] + (1 - global_param.ωL) * cp.Lᵉ
end


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
function update_π_cp!(
    cp::ConsumerGoodProducer
    )

    # total_freq = sum(map(machine -> machine.freq, cp.Ξ))
    cp.π = sum(map(machine -> (machine.freq * machine.A) / cp.n_machines, cp.Ξ))
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

    cp.mach_tb_repl = Vector{Machine}()

    for machine in cp.Ξ
        c_A = cp.w̄[end] / machine.A
        c_star = cp.w̄[end] / Aᵈ
        if c_A > c_star
            if p_star/(c_A - c_star) <= global_param.b || machine.age >= global_param.η
            # if p_star/(c_A - c_star) <= global_param.b
                push!(cp.mach_tb_repl, machine)
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
    Q = max(min(cp.π * cp.L, cp.π * cp.n_machines, cp.Qˢ), 0)
    # println("L ", cp.π * cp.L, " ", cp.n_machines, " ", cp.Qˢ)
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
        # c = wᵉ / cp.π
        c = cp.w̄[end] / cp.π
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

    # Check desired change in labor stock, also check for capital stock
    # as hiring more than this would not increase production.
    ΔLᵈ = min(cp.Qˢ, cp.n_machines * cp.π)/cp.π - cp.L

    if ΔLᵈ < -cp.L
        cp.ΔLᵈ = -cp.L
    else
        cp.ΔLᵈ = L * floor(ΔLᵈ / L)
    end
end


"""
Lets cp order machines from kp of choice.
"""
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


# """
# Computes profit of cp in previous time period
# """
# function compute_Π_cp!(
#     cp::ConsumerGoodProducer,
#     writeoffs=0.0::Float64
#     )
 
#     # Π = cp.curracc.S + cp.curracc.rev_dep - cp.curracc.TCL - cp.curracc.int_debt
#     # println(Π, " ", cp.curracc.S, " ", cp.curracc.TCL, " ", cp.curracc.TCI, " ", cp.curracc.int_debt)
#     Π = cp.curracc.S + cp.curracc.rev_dep - cp.curracc.TCL - cp.curracc.int_debt - writeoffs
#     push!(cp.Π, Π)
# end


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

    # # Check if complete order was received, otherwise partial replacement
    # if Iₜ == cp.Iᵈ
    #     # All investment received, remove all machines that were written off
    #     filter!(machine -> machine ∉ cp.mach_tb_repl, cp.Ξ)

    #     # Add new machine to capital stock
    #     push!(cp.Ξ, machine)
    # else
    #     # Partial investment received, remove part of machines.

    #     # TODO
    # end
    
    # All investment received, remove all machines that were written off
    filter!(machine -> machine ∉ cp.mach_tb_repl, cp.Ξ)

    # Add new machine to capital stock
    push!(cp.Ξ, machine)

    cp.curracc.TCI += Iₜ
end


# """
# Repays outstanding debts at start of period, shifts debt installments
# """
# function repay_debt_installments_cp!(
#     cp::ConsumerGoodProducer,
#     b=3::Float64
#     )

#     to_be_paid = cp.debt_installments[1]
#     cp.curracc.rep_debt += to_be_paid

#     cp.debt_installments[1] = cp.debt_installments[2]
#     cp.debt_installments[2] = cp.debt_installments[3]
#     cp.debt_installments[3] = 0.0
# end