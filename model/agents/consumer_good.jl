"""
Defines struct for consumer good producer
"""
mutable struct ConsumerGoodProducer <: AbstractAgent
    id :: Int                       # id
    type_good :: String             # type of good produced by producer
    age :: Int                      # firm age
    first_gen :: Bool               # shows if producer is in first generation 
    p :: Vector{Float64}            # hist prices
    c :: Vector{Float64}            # hist cost
    RD :: Vector{Float64}           # R&D spending
    D :: Vector{Float64}            # hist demand
    Dᵁ :: Float64                   # unsatisfied demand in last period
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
    n_machines :: Float64           # total freq of machines # TODO rename
    cu :: Float64                   # capital utilizataion
    employees:: Vector{Int}         # employees list
    L :: Float64                    # labor units
    Lᵉ:: Float64                    # exp labor force
    ΔLᵈ :: Float64                  # desired change in labor force
    w̄ :: Vector{Float64}            # wage level
    wᴼ :: Float64                   # offered wage
    wᴼ_max :: Float64               # maximum offered wage
    brochures :: Vector             # brochures from kp
    π :: Float64                    # hist productivity
    f :: Vector{Float64}            # hist market share
    μ :: Vector{Float64}            # markup rate
    Π :: Vector{Float64}            # hist profits
    NW_growth :: Float64            # growth rate of liquid assets
    cI :: Float64                   # internal funds for investments
    balance :: Balance              # balance sheet
    curracc :: FirmCurrentAccount   # current account
end


function initialize_cp(
    id::Int, 
    machines::Vector{Machine},  
    type_good::String,
    n_init_emp_cp::Int,
    μ1::Float64,
    ι::Float64;
    D=800.0::Float64,
    p=1.0 * (1 + μ1)::Float64,
    w=1.0::Float64,
    L=n_init_emp_cp * 100::Int,
    Lᵉ=n_init_emp_cp * 100::Int,
    N_goods=D*ι::Float64,
    n_consrgood=200::Int,
    f=2/n_consrgood,
    first_gen=true::Bool
    )

    balance = Balance(               
        0.0,                        # - N: inventory (in money units)
        0.0,                        # - K: capital
        0.0,                        # - NW: liquid assets
        0.0,                        # - debt: debt
        0.0                         # - EQ: equity
    )

    curracc = FirmCurrentAccount(0,0,0,0,0,0,0)

    cp = ConsumerGoodProducer(
        id,                         # global id
        type_good,                  # type of good produced by producer
        0,                          # firm age
        first_gen,
        [p],                        # p: hist prices
        [1.0],                      # c: hist cost
        [],                         # RD: hist R&D spending
        [D],                        # D: hist demand
        0.0,                        # Dᵁ unsatisfied demand in last period
        D,                          # Dᵉ exp demand
        Vector(),                   # hh_queue: vector containing ids of demanding households           
        ι * D,                      # Nᵈ: desired inventory
        N_goods,                    # N_goods: inventory in goods 
        [D * (1 + ι)],              # Q: hist production
        D * (1 + ι),                # Qᵉ: exp production
        0,                          # Qˢ: desired short-term production

        0,                          # Iᵈ: desired investments
        0,                          # EIᵈ: desired expansionary investments
        0,                          # RSᵈ: desired replacement investments
        [],                         # mach_tb_repl: list of to-be replaced machines
        0,                          # chosen_kp_id: id of chosen kp
        zeros(4),                   # debt installments

        machines,                   # Ξ: machines
        0,                          # n_machines: total number of machines
        0,                          # cu: capital utilization
        [],                         # employees: employees
        L,                          # L: labor units in company
        Lᵉ,                         # Lᵉ: exp labor force
        0.0,                        # ΔLᵈ: desired change in labor force
        [w],                        # w: wage level
        1.0,                        # wᴼ: offered wage
        1.0,                        # wᴼ_max: maximum offered wage
        [],                         # brochures from kp

        1.0,                        # π: hist productivity
        [f],                        # f: market share
        [μ1],                       # μ: hist markup
        [0.0],                      # Π: hist profits
        0.0,                        # NW_growth: growth rate of liquid assets
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
    update_μ_cp!(cp, global_param.μ1)

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
    update_Qᵉ_cp!(cp, global_param)

    # Update expected labor supply
    # update_Lᵉ_cp!(cp)
    # TODO: see if I still need this for something

    # Compute desired capital stock expansion
    # desired_n_machines = max(cp.Qᵉ / cp.π, cp.n_machines)
    desired_n_machines = max(cp.Qᵉ, cp.n_machines)
    # desired_n_machines = max(cp.Qˢ / cp.π, cp.n_machines)
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

            if cp.EIᵈ > max_add_debt
                # Decrease amount of expansionary investment.
                cp.Iᵈ = max_add_debt
                cp.EIᵈ = max_add_debt
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
            # cp.Qˢ = min((cp.L + cp.ΔLᵈ) * cp.π[end], cp.n_machines * cp.π[end])
        end
    end

    # Based on final production decisions, update max offered wage
    update_wᴼ_max_cp!(cp)
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

    Aᵈ = brochure[4]
    p_star = brochure[2]
    c_star = cp.w̄[end] / Aᵈ

    cp.mach_tb_repl = Vector{Machine}()

    for machine in cp.Ξ
        c_A = cp.w̄[end] / machine.A
        if (c_A != c_star && p_star/(c_A - c_star) <= global_param.b) || machine.age >= global_param.η
            push!(cp.mach_tb_repl, machine)
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
    if length(cp.brochures) == 0
        brochure = model[sample(all_kp)].brochure
        cp.chosen_kp_id = brochure[1]
        return brochure
    end

    # If multiple brochures, find brochure with lowest cost of production
    all_cop = map(broch -> broch[2] + b * cp.w̄[end] / broch[4], cp.brochures)
    brochure = cp.brochures[argmin(all_cop)]
    cp.chosen_kp_id = brochure[1]

    return brochure
end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(
    cp::ConsumerGoodProducer
    )

    # Compute total production amount
    # Q = max(min(cp.π * cp.L, cp.π * cp.n_machines, cp.Qˢ), 0)
    Q = max(min(cp.π * cp.L, cp.n_machines, cp.Qˢ), 0)
    # println("Q: $Q, πL: $(cp.π * cp.L), n_mach: $(cp.n_machines), Qˢ: $(cp.Qˢ), L: $(cp.L)")
    # @printf("L: %i, K: %i, Qˢ: %i\n", cp.L, cp.n_machines, cp.Qˢ)
    # @printf("πL: %i, πK: %i, Qˢ: %i\n", cp.π * cp.L, cp.π * cp.n_machines, cp.Qˢ)
    push!(cp.Q, Q)
    
    # Update rate of capital utilization
    if cp.n_machines > 0
        cp.cu = cp.Q[end] / cp.n_machines
    else
        cp.cu = 0
    end
    
    # Change inventory, will be amount households can buy from
    cp.N_goods += Q
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
        receive_order_kp!(model[cp.chosen_kp_id], cp.id, order)
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

    # Check if any orders received
    if length(cp.order_queue) == 0
        return
    end

    # Compute total demand for goods, regardless of if it can be fulfilled
    # Allows producers to produce more if production is insufficient
    total_D = sum(map(order -> order[2], cp.order_queue))
    push!(cp.D, total_D)

    total_unsat_demand = 0

    actual_S = 0.0

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
                total_unsat_demand += q - cp.N_goods
                cp.N_goods = 0.0
            else
                total_unsat_demand += q
            end

            tot_price = q * share_fulfilled * cp.p[end]
            actual_S += tot_price
            receive_order_hh!(model[hh_id], cp.id, tot_price, share_fulfilled)
        end
    end
    # println("Unsatisfied demand: ", total_unsat_demand)
    cp.Dᵁ = total_unsat_demand
    cp.curracc.S = actual_S
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


"""
Replaces cp, places cp in firm list of hh.
"""
function replace_bankrupt_cp!(
    bankrupt_bp::Vector{Int}, 
    bankrupt_lp::Vector{Int},
    bankrupt_kp::Vector{Int},
    all_hh::Vector{Int},
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    all_kp::Vector{Int},
    avg_cu::Float64,
    φ1::Float64,
    φ2::Float64,
    φ3::Float64,
    φ4::Float64,
    p̄::Float64,
    w̄::Float64,
    μ1::Float64,
    ι::Float64,
    model::ABM;
    n_machines_init=40::Int,
    n_init_emp_cp=10::Int
    )

    # Create vectors containing ids of non-bankrupt bp and lp
    nonbankrupt_bp = setdiff(all_bp, bankrupt_bp)
    nonbankrupt_lp = setdiff(all_lp, bankrupt_lp)

    avg_n_machines = mean(map(cp_id -> model[cp_id].n_machines, vcat(nonbankrupt_bp, nonbankrupt_lp)))
    avg_NW = mean(map(cp_id -> model[cp_id].balance.NW, vcat(nonbankrupt_bp, nonbankrupt_lp)))

    # Make weights for allocating cp to hh
    # Minimum is taken to avoid weird outcomes when all bp and lp went bankrupt
    weights_hh_bp = map(hh_id -> min(1, 1 / length(model[hh_id].bp)), all_hh)
    weights_hh_lp = map(hh_id -> min(1, 1 / length(model[hh_id].lp)), all_hh)

    # Generate a vector of all kp that did not go bankrupt
    # poss_kp = filter(kp_id -> kp_id ∉ bankrupt_kp, all_kp)

    for p_id in vcat(bankrupt_bp, bankrupt_lp)

        type_good="Basic"
        if p_id in bankrupt_lp
            type_good="Luxury"
        end

        # New cp receive a advanced type of machine, first select kp proportional
        # to their market share. cp can also select kp ids that went bankrupt in this 
        # period, as these producers have already been replaced with new companies
        weights_kp = map(kp_id -> model[kp_id].f[end], all_kp)
        kp_choice_id = sample(all_kp, Weights(weights_kp))
        A_choice = model[kp_choice_id].A[end]
        p_choice = model[kp_choice_id].p[end]

        # Sample what the size of the capital stock will be
        capital_coefficient = rand(Uniform(φ1, φ2))
        tot_freq_machines = capital_coefficient * avg_n_machines
        D = avg_cu * tot_freq_machines

        # Generate vector of machines
        machines = initialize_machine_stock(
                        tot_freq_machines, 
                        p_choice, 
                        A_choice, 
                        n_machines_init
                    )

        new_cp = initialize_cp(
                    p_id,
                    machines,
                    type_good,
                    0,
                    μ1,
                    ι;
                    D=D,
                    p=p̄,
                    w=w̄,
                    L=0,
                    Lᵉ=0,
                    N_goods=0.0,
                    first_gen=false
                )
        update_n_machines_cp!(new_cp)
        update_π_cp!(new_cp)

        # Sample what the stock of liquid assets will be
        NW_coefficient = rand(Uniform(φ3, φ4))
        NW_stock = NW_coefficient * avg_NW

        # Augment the balance with acquired NW and K
        new_cp.balance.NW = NW_stock
        new_cp.balance.K = tot_freq_machines * p_choice

        # Borrow funds for the machine and liquid assets
        # NW_to_borrow = NW_stock + tot_freq_machines * p_choice
        # borrow_funds_p!(new_cp, tot_freq_machines * p_choice + NW_to_borrow)
        # TODO: a kp firm has to receive revenue for this

        # new_cp.balance.NW = tot_freq_machines * p_choice
        # new_cp.balance.K = tot_freq_machines * p_choice
        # new_cp.balance.debt = tot_freq_machines * p_choice + NW_to_borrow

        add_agent!(new_cp, model)

        # Add new cp to subset of households, inversely proportional to amount of suppliers
        # they already have
        n_init_hh = 100
        # println("Len hh $(length(all_hh))")
        if p_id ∈ bankrupt_bp
            customers = sample(all_hh, Weights(weights_hh_bp), n_init_hh)
        else
            customers = sample(all_hh, Weights(weights_hh_lp), n_init_hh)
        end
    
        # Add cp to list of bp and lp, according to type
        for hh_id in customers
            if p_id ∈ bankrupt_bp
                push!(model[hh_id].bp, p_id)
            else
                push!(model[hh_id].lp, p_id)
            end
        end
    end
end


"""
UPDATING AND COMPUTING FUNCTIONS
"""

"""
Updates expected demand based Dᵉ
"""
function update_Dᵉ_cp!(
    cp::ConsumerGoodProducer,
    ωD::Float64
    )

    if length(cp.D) > 2
        # TODO: DESCRIBE IN MODEL 
        cp.Dᵉ = ωD * cp.Dᵉ + (1 - ωD) * (2 * cp.D[end] - cp.D[end-1])
    else
        cp.Dᵉ = ωD * cp.Dᵉ + (1 - ωD) * cp.D[end]
    end
end


"""
Updates desired short-term production Qˢ
"""
function update_Qˢ_cp!(
    cp::ConsumerGoodProducer
    )

    cp.Qˢ = max(cp.Dᵉ + cp.Nᵈ - cp.N_goods, 0)
end


"""
Updates expected long-term production Qᵉ
"""
function update_Qᵉ_cp!(
    cp::ConsumerGoodProducer,
    global_param::GlobalParam
    )

    # TODO: very large sensitivity to this parameter: check out!

    cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * (cp.Dᵉ + cp.Nᵈ)
    
    # if length(cp.D) > 2
    #     Qg = cp.Q[end] * (1 + (cp.D[end] - cp.D[end-1]) / cp.D[end-1])
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * Qg
    # else
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end]
    # end

    # cp.Qᵉ = cp.Dᵉ + cp.Dᵁ

    # if length(cp.Π) > 2 && cp.Π[end-1] != 0
    #     Qg = cp.Q[end] * (1 + (cp.Π[end] - cp.Π[end-1]) / cp.Π[end-1])
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * Qg
    # else
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end]
    # end

    # if cp.NW_growth != 0
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end] * (1 + cp.NW_growth)
    # else
    #     cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end]
    # end

    # cp.Qᵉ = global_param.ωQ * cp.Qᵉ + (1 - global_param.ωQ) * cp.Q[end] * (1 + cp.NW_growth)

    # cp.Qᵉ = cp.Q[end] + cp.Dᵁ

end


"""
Updates expected long-term labor supply Lᵉ
"""
function update_Lᵉ_cp!(
    cp::ConsumerGoodProducer,
    global_param
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

    cp.π = sum(map(machine -> (machine.freq * machine.A) / cp.n_machines, cp.Ξ))
end


"""
Computes the markup rate μ based on the market share f.
"""
function update_μ_cp!(
    cp::ConsumerGoodProducer, 
    μ1::Float64
    )

    # First check if market share is nonzero, else liquidate firm
    # if cp.f[end] <= 0.001
    #     # TODO: liquidate firm
    #     cp.f[end] = 0.001
    # end

    # # Compute market share
    # if length(cp.f) > 2
    #     cp.μ = cp.μ * (1 + υ * (cp.f[end] - cp.f[end-1])/cp.f[end-1])
    # else
    #     cp.μ = μ1
    # end
    # println(cp.μ)

    # println(cp.μ)
    # println(cp.Π)

    b = 0.1
    l = 4

    if length(cp.μ) > l && length(cp.Π) > l && cp.Π[end] != 0
        mean_μ = mean(cp.μ[end-l:end-1])
        # Δμ = (cp.μ[end] - cp.μ[end-1]) / cp.μ[end-1]
        Δμ = (cp.μ[end] - mean_μ) / mean_μ

        mean_Π = mean(cp.Π[end-l:end-1])
        # ΔΠ = (cp.Π[end] - cp.Π[end-1]) / cp.Π[end-1]
        ΔΠ = (cp.Π[end] - mean_Π) / mean_Π

        # println("$mean_μ, $mean_Π")
        # println("Δμ: $Δμ, $(sign(Δμ)), ΔΠ: $ΔΠ, $(sign(ΔΠ))")

        shock = rand(Uniform(0.0, b))

        new_μ = max(cp.μ[end] * (1 + sign(Δμ) * sign(ΔΠ) * shock), 0)
        push!(cp.μ, new_μ)

        # if ΔΠ > 0
        #     new_μ = cp.μ[end] * (1 + Δμ)
        #     push!(cp.μ, new_μ)
        # else
        #     new_μ = max(cp.μ[end] * (1 - Δμ), 0)
        #     push!(cp.μ, new_μ)
        # end

    elseif cp.Π[end] == 0
        push!(cp.μ, cp.μ[end] * (1 + rand(Uniform(-b, 0.0))))
    else
        push!(cp.μ, cp.μ[end] * (1 + rand(Uniform(-b, b))))
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
        # wᵉ = (cp.w̄[end] * (cp.L + cp.ΔLᵈ)) / (cp.L + cp.ΔLᵈ)
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

    p = (1 + cp.μ[end]) * cp.c[end]
    # println("p ", p, " c ", cp.c[end], " μ ", cp.μ)
    push!(cp.p, p)
end


"""
Computes the desired change in labor supply ΔLᵈ
    Check desired change in labor stock, also check for capital stock
    as hiring more than this would not increase production.
"""
function update_ΔLᵈ_cp!(
    cp::ConsumerGoodProducer
    )

    Lᵈ = min(cp.Qˢ / cp.π - cp.L, cp.n_machines / cp.π - cp.L)
    # println("Lᵈ: $(Lᵈ), L: $(cp.L), Qˢ/π: $(cp.Qˢ / cp.π) n_mach/π: $(cp.n_machines / cp.π)")
    cp.ΔLᵈ = max(Lᵈ, -cp.L)
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


"""
Updates maximum offered wage wᴼ_max
"""
function update_wᴼ_max_cp!(
    cp::ConsumerGoodProducer
    )
    # TODO: DESCRIBE IN MODEL
    cp.wᴼ_max = cp.π[end] * cp.p[end]
    # if cp.ΔLᵈ > 0
    #     # cp.wᴼ_max = (cp.Dᵉ * cp.p[end] - cp.w̄[end] * cp.L) / cp.ΔLᵈ
    #     cp.wᴼ_max = cp.p[end]
    # else
    #     cp.wᴼ_max = 0
    # end
end