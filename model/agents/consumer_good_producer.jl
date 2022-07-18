"""
Defines struct for consumer good producer
"""
@with_kw mutable struct ConsumerGoodProducer <: AbstractAgent
    id::Int                                   # id
    age::Int = 0                              # firm age
    t_next_update::Int                        # next update time
    
    # Price and cost data
    μ::Vector{Float64}                        # markup rate
    p::Vector{Float64} = fill(1+μ[end], 3)    # hist prices
    c::Vector{Float64} = ones(Float64, 3)     # hist cost
    true_c::Float64 = 0.0                     # true unit cost

    # Production and demand
    D::Vector{Float64}                        # hist demand
    Dᵁ::Vector{Float64} = zeros(3)            # unsatisfied demand in last period
    Dᵉ::Float64                               # exp demand
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
    mach_tb_retired::Vector{Machine} = []     # list of to-be retired machines (without replacement)
    chosen_kp_id::Int = 0                     # id of chosen kp
    debt_installments::Vector{Float64} = zeros(4) # installments of debt repayments

    Ξ::Vector{Machine}                        # machines
    n_machines::Float64 = 0                   # total freq of machines # TODO rename
    cu::Float64 = 0                           # capital utilizataion
    employees::Vector{Int} = []               # employees list
    L::Float64                                # labor units
    Lᵈ::Float64 = L                           # desired labor units
    ΔLᵈ::Float64 = 0.0                        # desired change in labor force
    w̄::Vector{Float64}                        # wage level
    wᴼ::Float64 = 1.0                         # offered wage
    wᴼ_max::Float64 = 1.0                     # maximum offered wage
    brochures::Vector = []                    # brochures from kp

    π_LP::Float64 = 1.0                       # labor productivity of total capital stock
    π_EE::Float64 = 1.0                       # productivity per energy unit of total capital stock
    π_EF::Float64 = 1.0                       # environmental friendlisness of total capital stock
    mean_skill::Float64 = 1.0                 # mean skill level of employees
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
    t_next_update::Int, 
    machines::Vector{Machine},  
    n_init_emp_cp::Int,
    μ::Float64,
    ι::Float64;
    D=1600.0::Float64,
    w=1.0::Float64,
    L=n_init_emp_cp * 100::Int,
    N_goods=D*ι::Float64,
    n_consrgood=200::Int,
    f=2/n_consrgood,
    )

    cp = ConsumerGoodProducer(
        id = id,
        t_next_update = t_next_update,
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

    cp.balance.NW = 1500
    cp.balance.EQ = 1500

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
    government::Government,
    ep,
    globalparam::GlobalParam,
    μ_avg::Float64,
    τˢ::Float64,
    t::Int,
    model::ABM
    )

    # Update amount of owned capital and desired inventories
    update_n_machines_cp!(cp, globalparam.freq_per_machine)
    update_Nᵈ_cp!(cp, globalparam.ι)

    # Compute expected demand
    update_Dᵉ_cp!(cp, globalparam.ω)

    # Compute desired short-term production
    update_Qˢ_cp!(cp)

    # Update average productivity
    update_π_cp!(cp)

    # Compute corresponding change in labor stock
    update_Lᵈ!(cp, globalparam.λ)

    # Update average wage w̄
    update_w̄_p!(cp, model)

    # Update markup μ
    update_μ_cp!(cp)

    # Update cost of production c
    compute_c_cp!(cp, ep.pₑ[t], government.τᴱ, government.τᶜ)

    # Compute price
    compute_p_cp!(cp, τˢ)
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
    government::Government, 
    all_kp::Vector{Int},
    globalparam::GlobalParam,
    ep,
    t::Int,
    model::ABM
    )

    # Set number of ordered machines to zero
    cp.n_mach_ordered_EI = 0
    cp.n_mach_ordered_RS = 0

    # Choose kp
    brochure = choose_producer_cp!(cp, government, globalparam.b, all_kp, ep, t, model)

    # Update LT production
    update_Qᵉ_cp!(cp, globalparam.ω, globalparam.ι)

    # Plan replacement investments
    plan_replacement_cp!(cp, government, globalparam, ep, brochure, t)

    # Plan expansion investments
    plan_expansion_cp!(cp, globalparam, brochure)

    # Determine total investments
    cp.Iᵈ = (cp.EIᵈ + cp.RSᵈ) / globalparam.update_period

    # See if enough funds available for investments and production, otherwise
    # change investments and production to match funding availability.

    check_funding_restrictions_cp!(cp, government, globalparam, ep.pₑ[t])

    # Send orders to kp
    order_machines_cp!(cp, model)
end


"""
Checks funding restructions based on expected revenue and expenses. If not enough
    funding available in firm, decrease desired production, hiring or investments.
"""
function check_funding_restrictions_cp!(
    cp::ConsumerGoodProducer,
    government::Government,
    globalparam::GlobalParam,
    pₑ::Float64
    )

    # Determine expected TCL and TCE
    TCLᵉ = cp.ΔLᵈ > 0 ? cp.w̄[end] * cp.L + cp.ΔLᵈ * cp.wᴼ : (cp.L + cp.ΔLᵈ) * cp.w̄[end]
    TCE = (1 + government.τᴱ) * pₑ * cp.Qˢ / cp.π_EE + government.τᶜ * cp.Qˢ * cp.π_EF

    # Determine how much additional debt can be made
    max_add_debt = max(globalparam.Λ * cp.D[end] * cp.p[end - 1] - cp.balance.debt, 0)

    # Check if cost of labor and investment can be financed from liquid assets
    NW_no_prod = (cp.balance.NW + cp.Dᵉ * cp.p[end] + cp.curracc.rev_dep 
                  - cp.debt_installments[1] - cp.balance.debt * globalparam.r)
    
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
    government::Government, 
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
    all_cop = Float64[]
    for brochure in cp.brochures
        p_mach = brochure[2]
        A_LP = brochure[4]
        A_EE = brochure[5]
        A_EF = brochure[6]

        c = p_mach + b * cop(cp.w̄[end], A_LP, government.τᴱ, ep.pₑ[t], A_EE, government.τᶜ, A_EF)

        push!(all_cop, c)
    end

    # Choose kp based on brochures
    brochure = sample(cp.brochures, Weights(1 ./ all_cop.^2))
    cp.chosen_kp_id = brochure[1]

    return brochure
end


"""
Plans replacement investment based on age machines and available new machines
"""
function plan_replacement_cp!(
    cp::ConsumerGoodProducer,
    government::Government,
    globalparam::GlobalParam,
    ep,
    brochure,
    t::Int
    )

    # Get price and cost of production of chosen kp
    p_star = brochure[2]
    c_star = cop(cp.w̄[end], brochure[4], government.τᴱ, ep.pₑ[t], brochure[5], government.τᶜ, brochure[6])

    cp.mach_tb_repl = []
    cp.mach_tb_retired = []

    # See if machine stock too large in order to decide if need to be replaced
    if cp.n_machines > cp.Qᵉ
        n_machines_too_many = cp.n_machines - cp.Qᵉ
    else
        n_machines_too_many = 0.0
    end

    # Sort machines by cost of production
    sort!(cp.Ξ , by = machine -> cp.w̄[end]/machine.A_LP + ep.pₑ[t]/machine.A_EE; rev=true)

    # Loop over machine stock, select which machines to replace
    for machine in cp.Ξ

        c_A = cop(cp.w̄[end], machine.A_LP, government.τᴱ, ep.pₑ[t], machine.A_EE, government.τᶜ, machine.A_EF)

        if machine.age >= globalparam.η
            # Machine has reached max age, decide if replaced or not
            if n_machines_too_many < machine.freq
                # No machines planned to be written off, replace old machine
                push!(cp.mach_tb_repl, machine)
            else
                # Do not replace machine
                push!(cp.mach_tb_retired, machine)
                n_machines_too_many -= machine.freq
            end

        elseif (c_A != c_star && p_star/(c_A - c_star) <= globalparam.b 
                && machine.age > globalparam.b)
            # New machine cheaper to operate, replace old machine
            push!(cp.mach_tb_repl, machine)
        end
    end

    # Sort to-be-replaces machines from least to most productive
    sort!(cp.mach_tb_repl, by=machine->machine.A_LP)

    # Update total amount of to-be-replaces machines
    cp.n_mach_ordered_RS = length(cp.mach_tb_repl)

    # Compute investment amount corresponding to replacement investments
    cp.RSᵈ = p_star * cp.n_mach_ordered_RS * globalparam.freq_per_machine
end


"""
Plans expansion investments based on expected production.
"""
function plan_expansion_cp!(
    cp::ConsumerGoodProducer,
    globalparam::GlobalParam,
    brochure
    )

    if cp.Qᵉ > cp.n_machines && cp.cu > 0.8
        cp.n_mach_ordered_EI = floor(Int64, (cp.Qᵉ - cp.n_machines) / globalparam.freq_per_machine)
        cp.EIᵈ = brochure[2] * cp.n_mach_ordered_EI * globalparam.freq_per_machine
    else
        cp.n_mach_ordered_EI = 0
        cp.EIᵈ = 0.0
    end
end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(
    cp::ConsumerGoodProducer,
    ep,
    globalparam::GlobalParam,
    t::Int
    )

    # If the cp does not need to use its complete capital stock, only use most productive 
    # machines
    n_machines_req = ceil(Int, cp.Qˢ / globalparam.freq_per_machine)
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
    Q = min(actual_π_LP * cp.L, cp.n_machines)
    shift_and_append!(cp.Q, Q)

    # Update energy use and carbon emissions from production
    update_EU_TCE_cp!(cp, actual_π_EE, ep.pₑ[t])
    update_emissions_cp!(cp, actual_em)
    
    # Update rate of capital utilization
    if cp.n_machines > 0
        cp.cu = Q / cp.n_machines
    else
        cp.cu = 0.0
    end
    
    # Change inventory, will be amount households can buy from
    cp.N_goods += Q
end


"""
    order_machines_cp!(cp::ConsumerGoodProducer, model::ABM)

Lets cp order machines from kp of choice.
"""
function order_machines_cp!(
    cp::ConsumerGoodProducer,
    model::ABM
    )

    total_n_machines = cp.n_mach_ordered_EI + cp.n_mach_ordered_RS

    # Send orders for machines to kp
    if total_n_machines > 0 && hascapacity(model[cp.chosen_kp_id])
        receive_order_kp!(model[cp.chosen_kp_id], cp.id, total_n_machines)
    end
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
function receive_machines_cp!(
    cp::ConsumerGoodProducer,
    ep, 
    new_machines::Vector{Machine},
    Iₜ::Float64,
    t::Int
    )

    cp.curracc.TCI += Iₜ

    # Replace old machines
    if length(cp.mach_tb_repl) > length(new_machines)
        # Not all to-be replaced machines were sent, only replace machines
        # that were delivered

        # Sort machines by cost of production, replace most expensive first
        sort!(cp.mach_tb_repl , by = machine -> cp.w̄[end]/machine.A_LP + ep.pₑ[t]/machine.A_EE; rev=true)
        filter!(machine -> machine ∉ cp.mach_tb_repl[1:length(new_machines)], cp.Ξ)
    else
        # All to-be replaced machines were sent, replace all machines and add
        # the additional machines as expansionary investment
        filter!(machine -> machine ∉ cp.mach_tb_repl, cp.Ξ)
    end

    cp.Ξ = vcat(cp.Ξ, new_machines)
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
    globalparam::GlobalParam,
    indexfund::IndexFund,
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
    capital_coefficients = rand(Uniform(globalparam.φ1, globalparam.φ2), n_bankrupt_cp)
    NW_coefficients = rand(Uniform(globalparam.φ3, globalparam.φ4), n_bankrupt_cp)

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
        all_n_machines[i] = floor(Int, capital_coefficients[i] * avg_n_machines / globalparam.freq_per_machine)
    end

    # Compute share of investments that can be paid from the investment fund
    req_NW = (avg_NW .* NW_coefficients) .+ (all_n_machines .* (kp_choice_ps .* globalparam.freq_per_machine))
    all_req_NW = sum(req_NW)
    frac_NW_if = decide_investments_if!(indexfund, all_req_NW, t)

    for (i,cp_id) in enumerate(bankrupt_cp)

        # Sample what the size of the capital stock will be
        D = macro_struct.cu[t] * all_n_machines[i] * globalparam.freq_per_machine

        # In the first period, the cp has no machines yet, these are delivered at the end
        # of the first period
        new_cp = initialize_cp(
                    cp_id,
                    t + 1,
                    Vector{Machine}(),
                    0,
                    macro_struct.μ_cp[t],     
                    globalparam.ι;
                    D=D,
                    w=macro_struct.w̄_avg[t],
                    L=0,
                    N_goods=0.0,
                    f=0.0
                )

        # Order machines at kp of choice
        new_cp.chosen_kp_id = kp_choice_ids[i]
        new_cp.n_mach_ordered_EI = all_n_machines[i]
        order_machines_cp!(new_cp, model)

        update_wᴼ_max_cp!(new_cp)

        # Augment the balance with acquired NW and K
        new_cp.balance.NW = req_NW[i]

        # Borrow remaining required funds for the machine, the other part of the 
        # funds come from the investment fund
        borrow_funds_p!(new_cp, (1 - frac_NW_if) * req_NW[i], globalparam.b)

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

    cp.Dᵉ = cp.age > 1 ? ω * cp.Dᵉ + (1 - ω) * (cp.D[end]) : cp.Dᵉ
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
    ω::Float64,
    ι::Float64
    )

    cp.Qᵉ = ω * cp.Qᵉ + (1 - ω) * ((1 + ι) * cp.Dᵉ)
end


"""
Updates capital stock n_machines
"""
function update_n_machines_cp!(
    cp::ConsumerGoodProducer,
    freq_per_machine::Int
    )

    # Retire old machines that are not replaced
    if length(cp.mach_tb_retired) > 0
        filter!(machine -> machine ∉ cp.mach_tb_retired, cp.Ξ)
    end

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
    cp.π_EF = length(cp.Ξ) > 0 ? sum(machine -> (machine.freq * machine.A_EF) / cp.n_machines, cp.Ξ) : 1.0
end


"""
Computes the markup rate μ based on the market share f.
"""
function update_μ_cp!(
    cp::ConsumerGoodProducer,
    )

    if cp.f[end] != 0.0 && cp.f[end-1] != 0.0
        new_μ = cp.μ[end] * min((1 + 0.04 * (cp.f[end] - cp.f[end-1]) / cp.f[end-1]), 1.04)
    else
        new_μ = cp.μ[end] * 0.95
    end
    shift_and_append!(cp.μ, new_μ)
end


"""
Compute production cost per unit c
"""
function compute_c_cp!(
    cp::ConsumerGoodProducer,
    pₑ::Float64,
    τᴱ::Float64,
    τᶜ::Float64
    )

    if cp.L + cp.ΔLᵈ > 0
        c = cop(cp.w̄[end], cp.π_LP, τᴱ, pₑ, cp.π_EE, τᶜ, cp.π_EF)
        shift_and_append!(cp.c, c)
    else
        shift_and_append!(cp.c, cp.c[end])
    end
end


"""
Computes price based on cost c and markup μ
"""
function compute_p_cp!(
    cp::ConsumerGoodProducer,
    τˢ::Float64
    )

    shift_and_append!(cp.p, (1 + τˢ)*((1 + cp.μ[end]) * max(cp.c[end], cp.true_c)))
end


"""
Computes the desired labor supply Lᵈ and the change in labor supply ΔLᵈ
    Check desired change in labor stock, also check for capital stock
    as hiring more than this would not increase production.
"""
function update_Lᵈ!(
    cp::ConsumerGoodProducer, 
    λ::Float64
    )

    cp.Lᵈ = λ * cp.L + (1 - λ) * min(cp.Qˢ / cp.π_LP, cp.n_machines / cp.π_LP)
    cp.ΔLᵈ = max(cp.Lᵈ - cp.L, -cp.L)
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
    # cp.wᴼ_max = (cp.π_LP * cp.p[end])
    cp.wᴼ_max = Inf
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