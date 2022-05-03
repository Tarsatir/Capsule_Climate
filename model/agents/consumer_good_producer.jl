"""
Defines struct for consumer good producer
"""
@with_kw mutable struct ConsumerGoodProducer <: AbstractAgent
    id::Int                                   # id
    age::Int = 0                              # firm age
    
    # Price and cost data
    μ::Vector{Float64}                        # markup rate
    p::Vector{Float64} = fill(1+μ[end], 3) # hist prices
    c::Vector{Float64} = ones(Float64, 3)     # hist cost

    # Production and demand
    D::Vector{Float64}                        # hist demand
    Dᵉ::Float64                               # exp demand
    Dᵁ::Float64 = 0.0                         # unsatisfied demand in last period
    order_queue::Vector = Vector()            # vector containing orders of demanding households
    Nᵈ::Float64                               # desired inventory
    N_goods::Float64                          # inventory in good units
    Q::Vector{Float64}                        # hist production
    Qᵉ::Float64                               # exp production
    Qˢ::Float64 = 0.0                         # desired short-term production
    EU::Float64 = 0.0                         # energy use in the last period

    # Investments
    Iᵈ::Float64 = 0                           # desired total investments
    EIᵈ::Float64 = 0                          # desired expansionary investments
    RSᵈ::Float64 = 0                          # desired replacement investments
    n_mach_ordered_EI::Int = 0                # number of machines ordered for expansion
    n_mach_ordered_RS::Int = 0                # number of machines ordered for replacement
    mach_tb_repl::Vector{Machine} = []        # list of to-be replaced machines
    chosen_kp_id::Int = 0                     # id of chosen kp
    debt_installments::Vector{Float64} = zeros(4) # installments of debt repayments

    Ξ::Vector{Machine}                        # machines
    n_machines::Float64 = 0                   # total freq of machines # TODO rename
    cu::Float64 = 0                           # capital utilizataion
    employees::Vector{Int} = []               # employees list
    L::Float64                                # labor units
    ΔLᵈ::Float64 = 0.0                        # desired change in labor force
    w̄::Vector{Float64}                        # wage level
    wᴼ::Float64 = 1.0                         # offered wage
    wᴼ_max::Float64 = 1.0                     # maximum offered wage
    brochures::Vector = []                    # brochures from kp

    π_LP::Float64 = 1.0                       # labor productivity of total capital stock
    π_EE::Float64 = 1.0                       # productivity per energy unit of total capital stock
    f::Vector{Float64}                        # hist market share
    Π::Vector{Float64} = zeros(Float64, 3)    # hist profits
    Πᵀ::Vector{Float64} = zeros(Float64, 3)   # Historical profits after tax
    NW_growth::Float64 = 0.0                  # growth rate of liquid assets
    cI::Float64 = 0.0                         # internal funds for investments
    balance::Balance = Balance()              # balance sheet
    curracc::FirmCurrentAccount = FirmCurrentAccount() # current account

    emissions::Float64 = 0.0                  # carbon emissions in last period
end


function initialize_cp(
    id::Int, 
    machines::Vector{Machine},  
    # type_good::String,
    n_init_emp_cp::Int,
    μ::Float64,
    ι::Float64;
    D=800.0::Float64,
    w=1.0::Float64,
    L=n_init_emp_cp * 100::Int,
    N_goods=3*D*ι::Float64,
    n_consrgood=200::Int,
    f=2/n_consrgood,
    )


    cp = ConsumerGoodProducer(
        id=id,
        μ = fill(μ, 3),
        D = fill(D, 3),
        Dᵉ = D,  
        Nᵈ = ι * D,                
        N_goods = N_goods,          
        Q = fill(D * (1 + ι), 3),   
        Qᵉ = D * (1 + ι),          
        Ξ = machines,                 
        L = L,
        w̄ = fill(w, 3),
        f = fill(f, 3)
    )

    cp.balance.NW = 500
    cp.balance.EQ = 500

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
    t::Int,
    model::ABM
    )

    # Update amount of owned capital and desired inventories
    update_n_machines_cp!(cp, global_param.freq_per_machine)
    update_Nᵈ_cp!(cp, global_param.ι)

    # Compute expected demand
    update_Dᵉ_cp!(cp, global_param.ω)

    # Compute desired short-term production
    update_Qˢ_cp!(cp)

    # Update average productivity
    update_π_cp!(cp)

    # Compute corresponding change in labor stock
    update_ΔLᵈ_cp!(cp)

    # Update average wage w̄
    update_w̄_p!(cp, model)

    # Update markup μ
    update_μ_cp!(cp, t, global_param.μ1)

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
    all_kp::Vector{Int},
    global_param::GlobalParam,
    ep,
    t::Int,
    model::ABM
    )

    # Set number of ordered machines to zero
    cp.n_mach_ordered_EI = 0
    cp.n_mach_ordered_RS = 0

    # Choose kp
    brochure = choose_producer_cp!(cp, global_param.b, all_kp, ep, t, model)

    # Update LT production
    update_Qᵉ_cp!(cp, global_param.ω)

    # Plan replacement investments
    plan_replacement_cp!(cp, global_param, ep, brochure, t)

    # Plan expansion investments
    plan_expansion_cp!(cp, global_param, brochure)

    # Determine total investments
    cp.Iᵈ = cp.EIᵈ + cp.RSᵈ

    # See if enough funds available for investments and production, otherwise
    # change investments and production to match funding availability.
    check_funding_restrictions_cp!(cp, global_param, ep.pₑ[t])

    # Send orders to kp
    order_machines_cp!(cp, model)
end


"""
Checks funding restructions based on expected revenue and expenses. If not enough
    funding available in firm, decrease desired production, hiring or investments.
"""
function check_funding_restrictions_cp!(
    cp::ConsumerGoodProducer,
    global_param::GlobalParam,
    pₑ::Float64
    )

    # Determine expected TCL and TCE
    TCLᵉ = cp.ΔLᵈ > 0 ? cp.w̄[end] * cp.L + cp.ΔLᵈ * cp.wᴼ : (cp.L + cp.ΔLᵈ) * cp.w̄[end]
    TCE = pₑ * cp.Qˢ / cp.π_EE

    # Determine how much additional debt can be made
    max_add_debt = max(global_param.Λ * cp.D[end] * cp.p[end - 1] - cp.balance.debt, 0)
    # max_add_debt = max(global_param.Λ * cp.Dᵉ * cp.p[end - 1] - cp.balance.debt, 0)

    # Check if cost of labor and investment can be financed from liquid assets
    NW_no_prod = (cp.balance.NW + cp.Dᵉ * cp.p[end] + cp.curracc.rev_dep 
                  - cp.debt_installments[1] - cp.balance.debt * global_param.r)
    
    if NW_no_prod > TCLᵉ + TCE + cp.Iᵈ
        # All cost of costs can be paid from liquid assets. No additional debt needed.
        cp.cI = cp.Iᵈ

    elseif NW_no_prod > TCLᵉ + TCE && NW_no_prod - TCLᵉ - TCE < cp.Iᵈ
        # Cost of labor can be paid from liquid assets, investment has to be partially 
        # funded from debt.
        cp.cI = NW_no_prod - TCLᵉ - TCE
        req_debt = cp.Iᵈ - (NW_no_prod - TCLᵉ - TCE)

        # Check if investment can be financed from additional debt, otherwise decrease investments
        if req_debt > max_add_debt
            if cp.EIᵈ > cp.cI + max_add_debt

                # Decrease amount of expansionary investment.
                poss_EI = cp.cI + max_add_debt
                cp.n_mach_ordered_EI = floor(Int, cp.n_mach_ordered_EI * (poss_EI / cp.EIᵈ))
                cp.n_mach_ordered_RS = 0
                cp.mach_tb_repl = []

            else

                # Full expansion is possible, decrease amount of replacement investments
                poss_RS = cp.cI + max_add_debt - cp.EIᵈ
                cp.n_mach_ordered_RS = floor(Int, cp.n_mach_ordered_RS * (poss_RS / cp.RSᵈ))

                if cp.n_mach_ordered_RS > 0
                    cp.mach_tb_repl = cp.mach_tb_repl[1:cp.n_mach_ordered_RS]
                else
                    cp.mach_tb_repl = []
                end

            end
        end

    else
        # Cost of labor exceeds liquid assets. Check if enough additional debt available.
        # All investment cancelled.
        cp.cI = 0.0
        cp.n_mach_ordered_EI = 0
        cp.n_mach_ordered_RS = 0
        cp.mach_tb_repl = Vector{Machine}()

        if NW_no_prod + max_add_debt < TCLᵉ + TCE
            # Cost of labor exceeds expected liquid assets plus max additional debt. 
            # Decrease production quantity.

            poss_prod = (NW_no_prod + max_add_debt) / (cp.w̄[end] / cp.π_LP + pₑ / cp.π_EE)
            poss_L = poss_prod / cp.π_LP
            cp.ΔLᵈ = poss_L - cp.L
            # println("Yeet: $(cp.ΔLᵈ)")
        end
    end

    # Based on final production decisions, update max offered wage
    update_wᴼ_max_cp!(cp)
end


"""
Lets cp make decision for kp out of available kp in brochures.
"""
function choose_producer_cp!(
    cp::ConsumerGoodProducer, 
    b::Int, 
    all_kp::Vector{Int},
    ep,
    t::Int,
    model::ABM
    )

    # In case of no brochures, pick a random kp
    if length(cp.brochures) == 0
        brochure = model[sample(all_kp)].brochure
        cp.chosen_kp_id = brochure[1]
        return brochure
    end

    # If multiple brochures, find brochure with lowest cost of production
    all_cop = map(broch -> broch[2] + b * (cp.w̄[end] / broch[4] + ep.pₑ[t] / broch[5]), cp.brochures)
    brochure = cp.brochures[argmin(all_cop)]
    cp.chosen_kp_id = brochure[1]

    return brochure
end


"""
Plans replacement investment based on age machines and available new machines
"""
function plan_replacement_cp!(
    cp::ConsumerGoodProducer,
    global_param::GlobalParam,
    ep,
    brochure,
    t::Int
    )

    # Aᵈ_LP = brochure[4]
    p_star = brochure[2]
    c_star = cp.w̄[end] / brochure[4] + ep.pₑ[t] / brochure[5]

    # Loop over machine stock, select which machines to replaceF
    cp.mach_tb_repl = []
    for machine in cp.Ξ

        c_A = cp.w̄[end] / machine.A_LP + ep.pₑ[t] / brochure[5]

        if ((c_A != c_star && p_star/(c_A - c_star) <= global_param.b) 
            || machine.age >= global_param.η)
            push!(cp.mach_tb_repl, machine)
        end
    end

    # Sort to-be-replaces machines from least to most productive
    sort!(cp.mach_tb_repl, by=machine->machine.A_LP)

    # Update total amount of to-be-replaces machines
    cp.n_mach_ordered_RS = length(cp.mach_tb_repl)

    # Compute investment amount corresponding to replacement investments
    cp.RSᵈ = p_star * cp.n_mach_ordered_RS * global_param.freq_per_machine
end


"""
Plans expansion investments based on expected production.
"""
function plan_expansion_cp!(
    cp::ConsumerGoodProducer,
    global_param::GlobalParam,
    brochure
    )

    if cp.Qᵉ > cp.n_machines
        cp.n_mach_ordered_EI = floor((cp.Qᵉ - cp.n_machines) / global_param.freq_per_machine)
        cp.EIᵈ = brochure[2] * cp.n_mach_ordered_EI * global_param.freq_per_machine
    else
        cp.EIᵈ = 0
    end
end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(
    cp::ConsumerGoodProducer,
    ep,
    global_param::GlobalParam,
    t::Int
    )

    # If the cp does not need to use its complete capital stock, only use most productive 
    # machines
    n_machines_req = ceil(Int, cp.Qˢ / global_param.freq_per_machine)
    if n_machines_req < length(cp.Ξ)
        # Compute number of machines needed (machines already ordered on productivity, 
        # least to most productive)
        req_machines = cp.Ξ[end-n_machines_req:end]
        actual_π_LP = mean(machine -> machine.A_LP, req_machines)
        actual_π_EE = mean(machine -> machine.A_EE, req_machines)
        actual_em = mean(machine -> machine.A_EF, req_machines)
    else
        actual_π_LP = cp.π_LP[end]
        actual_π_EE = cp.π_EE[end]
        actual_em = length(cp.Ξ) > 0 ? mean(machine -> machine.A_EF, cp.Ξ) : 0.0
    end

    # Compute total production amount
    Q = max(min(actual_π_LP * cp.L, cp.n_machines, cp.Qˢ), 0)
    shift_and_append!(cp.Q, Q)

    # Update energy use and carbon emissions from production
    update_EU_TCE_cp!(cp, actual_π_EE, ep.pₑ[t])
    update_emissions_cp!(cp, actual_em)
    
    # Update rate of capital utilization
    if cp.n_machines > 0
        cp.cu = Q / cp.n_machines
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

    total_n_machines = cp.n_mach_ordered_EI + cp.n_mach_ordered_RS

    # Send orders for machines to kp
    for _ in 1:total_n_machines
        receive_order_kp!(model[cp.chosen_kp_id], cp.id)
    end
end


# """
# Decides if enough inventory to send orders hh, 
#     if yes: transact, if no: send back order unfulfilled or partially fulfilled
# """
# function send_ordered_goods_all_cp!(
#     all_cp::Vector{Int},
#     t::Int,
#     model::ABM,
#     to
#     )

#     # Keep track of total demand, sales and unsatisfied demand
#     total_D = 0.0
#     actual_S = 0.0
#     total_unsat_demand = 0.0
#     share_fulfilled::Float64 = 0.0

#     for cp_id in all_cp
#         # Check if any orders received
#         if model[cp_id].age == 1 && t != 1
#             return nothing
#         elseif length(model[cp_id].order_queue) == 0
#             shift_and_append!(model[cp_id].D, 0.0)
#             model[cp_id].Dᵁ = 0.0
#             model[cp_id].curracc.S = 0.0
#             return nothing
#         end

#         # Keep track of total demand, sales and unsatisfied demand
#         total_D = 0.0
#         actual_S = 0.0
#         total_unsat_demand = 0.0
#         share_fulfilled = 0.0
        
#         # Loop over orders in queue, add realized sales S
#         for (hh_id, q) in model[cp_id].order_queue

#             total_D += q

#             # Check if enough inventory available
#             if model[cp_id].N_goods > q
#                 share_fulfilled = 1.0
#                 model[cp_id].N_goods -= q
#             elseif model[cp_id].N_goods > 0 && model[cp_id].N_goods < q
#                 share_fulfilled = ceil(model[cp_id].N_goods / q, digits=2)
#                 total_unsat_demand += q - model[cp_id].N_goods
#                 model[cp_id].N_goods = 0
#             else
#                 share_fulfilled = 0.0
#                 total_unsat_demand += q
#             end

#             if share_fulfilled != 1.0
#                 println(share_fulfilled)
#             end

#             tot_price = q * share_fulfilled * model[cp_id].p[end]
#             actual_S += tot_price
#             @timeit to "receive order" receive_ordered_goods_hh!(model[hh_id], model[cp_id].id, tot_price, share_fulfilled)
#         end

#         shift_and_append!(model[cp_id].D, total_D)
#         model[cp_id].Dᵁ = total_unsat_demand
#         model[cp_id].curracc.S = actual_S

#         reset_queue_cp!(model[cp_id])
#     end
# end


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
function receive_machines_cp!(
    cp::ConsumerGoodProducer, 
    new_machines::Vector{Machine}, 
    Iₜ::Float64,
    )

    # Add new machines to capital stock, add investment expenditure
    cp.Ξ = vcat(cp.Ξ, new_machines)
    cp.curracc.TCI += Iₜ

    # Check if complete order was received, otherwise partial replacement
    if cp.n_mach_ordered_EI == 0 && cp.n_mach_ordered_RS > length(new_machines)

        # All machines are meant for replacement
        cp.mach_tb_repl = cp.mach_tb_repl[1:length(new_machines)]

    elseif length(new_machines) < cp.n_mach_ordered_EI + cp.n_mach_ordered_RS

        # Decide how many old machines to discard
        n_old_replace = max(length(new_machines) - cp.n_mach_ordered_EI, 0)
        if n_old_replace > 0
            cp.mach_tb_repl = cp.mach_tb_repl[1:n_old_replace+1]
        else
            cp.mach_tb_repl = []
        end
    end

    filter!(machine -> machine ∉ cp.mach_tb_repl, cp.Ξ)
end


"""
Replaces cp, places cp in firm list of hh.
"""
function replace_bankrupt_cp!(
    bankrupt_cp::Vector{Int},
    bankrupt_kp::Vector{Int},
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    global_param::GlobalParam,
    indexfund_struct::IndexFund,
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    # Create vectors containing ids of non-bankrupt bp, lp and kp
    nonbankrupt_cp = setdiff(all_cp, bankrupt_cp)
    nonbankrupt_kp = setdiff(all_kp, bankrupt_kp)

    avg_n_machines = mean(cp_id -> model[cp_id].n_machines, nonbankrupt_cp)
    avg_NW = mean(cp_id -> model[cp_id].balance.NW, nonbankrupt_cp)

    # Make weights for allocating cp to hh
    # Minimum is taken to avoid weird outcomes when all bp and lp went bankrupt
    weights_hh_cp = map(hh_id -> min(1, 1 / length(model[hh_id].cp)), all_hh)
    weights_kp = map(kp_id -> model[kp_id].f[end], nonbankrupt_kp)

    n_bankrupt_cp = length(bankrupt_cp)

    # Sample all NW coefficients and capital coefficients
    capital_coefficients = rand(Uniform(global_param.φ1, global_param.φ2), n_bankrupt_cp)
    NW_coefficients = rand(Uniform(global_param.φ3, global_param.φ4), n_bankrupt_cp)

    # New cp receive a advanced type of machine, first select kp proportional
    # to their market share. cp can also select kp ids that went bankrupt in this 
    # period, as these producers have already been replaced with new companies
    kp_choice_ids = zeros(Int, n_bankrupt_cp)
    kp_choice_ps = zeros(Float64, n_bankrupt_cp)
    all_n_machines = zeros(Int, n_bankrupt_cp)
    for i in 1:n_bankrupt_cp
        # Decide from which kp to buy
        kp_choice_ids[i] = sample(all_kp, Weights(weights_kp))
        kp_choice_ps[i] = model[kp_choice_ids[i]].p[end]

        # Compute the number of machines each cp will buy
        all_n_machines[i] = floor(Int, capital_coefficients[i] * avg_n_machines / global_param.freq_per_machine)
    end

    # Compute share of investments that can be paid from the investment fund
    req_NW = (avg_NW .* NW_coefficients) .+ (all_n_machines .* (kp_choice_ps .* global_param.freq_per_machine))
    all_req_NW = sum(req_NW)
    frac_NW_if = decide_investments_if!(indexfund_struct, all_req_NW, t)

    for (i,cp_id) in enumerate(bankrupt_cp)

        # Sample what the size of the capital stock will be
        D = macro_struct.cu[t] * all_n_machines[i] * global_param.freq_per_machine
        # println("D: $D, cu: $(macro_struct.cu[t]), am: $(all_n_machines[i])")

        # In the first period, the cp has no machines yet, these are delivered at the end
        # of the first period
        new_cp = initialize_cp(
                    cp_id,
                    Vector{Machine}(),
                    0,
                    macro_struct.μ_cp[t],     
                    global_param.ι;
                    D=D,
                    w=macro_struct.w̄_avg[t],
                    L=0,
                    # Lᵉ=0,
                    N_goods=0.0,
                    f=0.0
                )

        # Order machines at kp of choice
        new_cp.chosen_kp_id = kp_choice_ids[i]
        new_cp.n_mach_ordered_EI = all_n_machines[i]
        order_machines_cp!(new_cp, model)

        # Augment the balance with acquired NW and K
        new_cp.balance.NW = req_NW[i]

        # Borrow remaining required funds for the machine, the other part of the 
        # funds come from the investment fund
        borrow_funds_p!(new_cp, (1 - frac_NW_if) * req_NW[i], global_param.b)

        add_agent!(new_cp, model)

        # Add new cp to subset of households, inversely proportional to amount of suppliers
        # they already have
        n_init_hh = 100

        customers = sample(all_hh, Weights(weights_hh_cp), n_init_hh)
    
        # Add cp to list of hh
        for hh_id ∈ customers
            push!(model[hh_id].cp, cp_id)
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
    ω::Float64
    )

    if cp.age > 2
        # TODO: DESCRIBE IN MODEL 
        cp.Dᵉ = ω * cp.Dᵉ + (1 - ω) * (2 * cp.D[end] - cp.D[end-1])
    elseif cp.age == 2
        cp.Dᵉ = ω * cp.Dᵉ + (1 - ω) * cp.D[end]
    end
end


"""
Updates desired short-term production Qˢ
"""
function update_Qˢ_cp!(
    cp::ConsumerGoodProducer
    )

    cp.Qˢ = max(cp.Dᵉ + cp.Nᵈ - cp.N_goods, 0.0)
end


"""
Updates expected long-term production Qᵉ
"""
function update_Qᵉ_cp!(
    cp::ConsumerGoodProducer,
    ω::Float64
    )

    cp.Qᵉ = ω * cp.Qᵉ + (1 - ω) * (cp.Dᵉ + cp.Nᵈ)
end


"""
Updates capital stock n_machines
"""
function update_n_machines_cp!(
    cp::ConsumerGoodProducer,
    freq_per_machine::Int
    )

    cp.n_machines = length(cp.Ξ) * freq_per_machine
end


"""
Updates weighted producivity of machine stock π_LP
"""
function update_π_cp!(
    cp::ConsumerGoodProducer
    )

    cp.π_LP = length(cp.Ξ) > 0 ? sum(machine -> (machine.freq * machine.A_LP) / cp.n_machines, cp.Ξ) : 1.0
    cp.π_EE = length(cp.Ξ) > 0 ? sum(machine -> (machine.freq * machine.A_EE) / cp.n_machines, cp.Ξ) : 1.0
end


"""
Computes the markup rate μ based on the market share f.
"""
function update_μ_cp!(
    cp::ConsumerGoodProducer,
    t::Int, 
    μ1::Float64
    )

    b = 0.1
    l = 2
    new_μ = cp.μ[end]

    # TODO describe Calvo Pricing

    # if rand() < 1/3

    if cp.age > l && cp.Π[end] != 0

        mean_μ = mean(cp.μ)
        Δμ = (cp.μ[end] - mean_μ) / mean_μ

        mean_Π = mean(cp.Π)
        ΔΠ = (cp.Π[end] - mean_Π) / mean_Π

        # TODO: CHANGE TO TAXED PROFITS??

        shock = rand(Uniform(0.0, b))
        new_μ = max(cp.μ[end] * (1 + sign(Δμ) * sign(ΔΠ) * shock), 0)

    elseif cp.Π[end] == 0
        new_μ = cp.μ[end] * (1 + rand(Uniform(-b, 0.0)))
    else
        new_μ = cp.μ[end] * (1 + rand(Uniform(-b, b)))
    end

    shift_and_append!(cp.μ, new_μ)
    # else
        # shift_and_append!(cp.μ, cp.μ[end])
    # end
end


"""
Compute production cost per unit c
"""
function compute_c_cp!(
    cp::ConsumerGoodProducer,
    )

    if cp.L + cp.ΔLᵈ > 0
        c = cp.w̄[end] / cp.π_LP
        shift_and_append!(cp.c, c)
    else
        shift_and_append!(cp.c, cp.c[end])
    end
end


"""
Computes price based on cost c and markup μ
"""
function compute_p_cp!(
    cp::ConsumerGoodProducer
    )

    shift_and_append!(cp.p, (1 + cp.μ[end]) * cp.c[end])
end


"""
Computes the desired change in labor supply ΔLᵈ
    Check desired change in labor stock, also check for capital stock
    as hiring more than this would not increase production.
"""
function update_ΔLᵈ_cp!(
    cp::ConsumerGoodProducer
    )

    Lᵈ = min(cp.Qˢ / cp.π_LP - cp.L, cp.n_machines / cp.π_LP - cp.L)
    cp.ΔLᵈ = max(Lᵈ, -cp.L)
end


"""
Updates desired inventory to be a share of the previous demand
"""
function update_Nᵈ_cp!(
    cp::ConsumerGoodProducer,
    ι::Float64
    )

    if cp.age > 1
        cp.Nᵈ = ι * cp.D[end]
    end
end


"""
Updates maximum offered wage wᴼ_max
"""
function update_wᴼ_max_cp!(
    cp::ConsumerGoodProducer
    )
    # TODO: DESCRIBE IN MODEL
    cp.wᴼ_max = cp.π_LP * cp.p[end]
    # if cp.ΔLᵈ > 0
    #     # cp.wᴼ_max = (cp.Dᵉ * cp.p[end] - cp.w̄[end] * cp.L) / cp.ΔLᵈ
    #     cp.wᴼ_max = cp.p[end]
    # else
    #     cp.wᴼ_max = 0
    # end
end


"""
Updates energy use for production
"""
function update_EU_TCE_cp!(
    cp::ConsumerGoodProducer,
    actual_π_EE::Float64,
    pₑ::Float64
    )

    cp.EU = length(cp.Ξ) > 0 ? cp.Q[end] / actual_π_EE : 0.0
    cp.curracc.TCE = pₑ * cp.EU
end


"""
Updates carbon emissions during production
"""
function update_emissions_cp!(
    cp::ConsumerGoodProducer, 
    actual_em::Float64
    )

    cp.emissions = actual_em * cp.EU
end