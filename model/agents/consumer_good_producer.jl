"""
Defines struct for consumer good producer
"""
@with_kw mutable struct ConsumerGoodProducer <: AbstractAgent

    id::Int64                                 # global agent id
    cp_i::Int64                               # cp index
    age::Int64 = 0                            # firm age
    # t_next_update::Int64                      # next update time
    
    # Price and cost data
    μ::Vector{Float64}                        # markup rate
    p::Vector{Float64} = fill(1+μ[end], 3)    # hist prices
    c::Vector{Float64} = ones(Float64, 3)     # hist cost
    true_c::Float64 = 0.                      # true unit cost

    # Production and demand
    D::Vector{Float64}                        # hist demand
    Dᵁ::Vector{Float64} = zeros(Float64, 3)   # unsatisfied demand in last period
    Dᵉ::Float64                               # exp demand
    order_queue::Vector = Vector()            # vector containing orders of demanding households
    Nᵈ::Float64                               # desired inventory
    N_goods::Float64                          # inventory in good units
    Q::Vector{Float64}                        # hist production
    Qᵉ::Float64                               # exp production
    Qˢ::Float64 = 0.                          # desired short-term production
    EU::Float64 = 0.                          # energy use in the last period

    # Investments
    possible_I::Float64 = 0.                  # Available funds for investments
    Iᵈ::Float64 = 0.                          # desired total investments
    EIᵈ::Float64 = 0.                         # desired expansionary investments
    RSᵈ::Float64 = 0.                         # desired replacement investments
    n_mach_desired_EI::Int64 = 0              # number of machines desired for expansion
    n_mach_ordered_EI::Int64 = 0              # number of machines ordered for expansion
    n_mach_desired_RS::Int64 = 0              # number of machines desired for replacement
    n_mach_ordered_RS::Int64 = 0              # number of machines ordered for replacement
    n_mach_desired_total::Int64 = 0           # number of total machines desired
    mach_tb_repl::Vector{Machine} = []        # list of to-be replaced machines
    mach_tb_retired::Vector{Machine} = []     # list of to-be retired machines (without replacement)
    kp_ids::Vector{Int64} = zeros(Int64, 0)   # ids of known kp
    mach_desired_per_kp::SortedDict{Int64, Tuple{Int64, Int64}} = SortedDict()
    debt_installments::Vector{Float64}        # installments of debt repayments

    Ξ::Vector{Machine}                        # machines
    n_machines::Float64 = 0.                  # total freq of machines # TODO rename
    cu::Float64 = 0.                          # capital utilizataion
    employees::Vector{Int64} = []               # employees list
    L::Float64                                # labor units
    Lᵈ::Float64 = L                           # desired labor units
    ΔLᵈ::Float64 = 0.0                        # desired change in labor force
    w̄::Vector{Float64}                        # wage level
    wᴼ::Float64 = 1.0                         # offered wage
    wᴼ_max::Float64 = 1.0                     # maximum offered wage
    # brochures::Vector = []                    # brochures from kp

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
    id::Int64,
    cp_i::Int64,
    # t_next_update::Int64, 
    machines::Vector{Machine},
    model::ABM;
    D::Float64 = 1600.,
    w::Float64 = 1.,
    f::Float64 = 1 / model.i_param.n_cp,
)

    cp = ConsumerGoodProducer(
        id = id,
        cp_i = cp_i,
        # t_next_update = t_next_update,
        μ = fill(model.g_param.μ1, 3),
        D = fill(D, 3),
        Dᵉ = D,  
        Nᵈ = model.g_param.ι * D,                
        N_goods = D * model.g_param.ι,          
        Q = fill(D * (1 + model.g_param.ι), 3),   
        Qᵉ = D * (1 + model.g_param.ι),
        debt_installments = zeros(Float64, model.g_param.b+1),          
        Ξ = machines,                 
        L = 0,
        w̄ = fill(w, 3),
        f = fill(f, 3)
    )

    cp.balance.NW = 1500.
    cp.balance.EQ = 1500.

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
    τˢ::Float64,
    n_hh::Int64,
    n_cp::Int64,
    t::Int64,
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

    # Update cost of production c
    compute_c_cp!(cp, ep.p_ep[t], government.τᴱ, government.τᶜ)

    # if checkupdateprice(cp.id, t, n_hh, n_cp, globalparam.p_rigid_time)
    if rand() < 1 / globalparam.p_rigid_time

        # Update markup μ
        update_μ_p!(cp, globalparam.ϵ_μ, t)

        # Compute price
        compute_p_cp!(cp, τˢ)
    end
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
    all_kp::Vector{Int64},
    globalparam::GlobalParam,
    ep::AbstractAgent,
    t::Int64,
    model::ABM
)

    # Rank producers
    rank_producers_cp!(cp, government, globalparam.b, all_kp, ep, t, model)

    # Update LT production
    update_Qᵉ_cp!(cp, globalparam.ω, globalparam.ι)

    # Plan replacement investments
    plan_replacement_cp!(cp, government, globalparam, ep, t, model)

    # Plan expansion investments
    plan_expansion_cp!(cp, globalparam, model)

    # Determine total investments
    cp.Iᵈ = cp.EIᵈ + cp.RSᵈ

    # See if enough funds available for investments and production, otherwise
    # change investments and production to match funding availability.
    check_funding_restrictions_cp!(cp, government, globalparam, ep.p_ep[t])

    # Send orders to kp
    order_machines_cp!(cp, globalparam.freq_per_machine, model)
end


# function plan_investment_cp!(
#     cp::ConsumerGoodProducer,
#     government::Government,
#     globalparam::GlobalParam,
#     ep::AbstractAgent,
#     all_kp::Vector{Int64},
#     t::Int64,
#     model::ABM
# )

#     # Rank machines based on cop
#     rank_machines_cp!(cp)

#     # Rank producers
#     rank_producers_cp!(cp, government, globalparam.b, all_kp, ep, t, model)

#     # Determine how many machines are desired per known capital good producer
#     # cp.mach_desired_per_kp = SortedDict(kp_id => (plan_replacement_cp!(), plan_expansion_cp!()) for kp_id in cp.kp_ids)
# end


"""
Checks funding restructions based on expected revenue and expenses. If not enough
    funding available in firm, decrease desired production, hiring or investments.
"""
function check_funding_restrictions_cp!(
    cp::ConsumerGoodProducer,
    government::Government,
    globalparam::GlobalParam,
    p_ep::Float64
    )

    # Determine expected TCL and TCE
    TCLᵉ = (cp.L + cp.ΔLᵈ) * cp.w̄[end]
    TCE = cp.Qˢ * ((p_ep + government.τᴱ) / cp.π_EE + government.τᶜ * cp.π_EF)

    # Determine how much additional debt can be made
    max_add_debt = max(globalparam.Λ * cp.D[end] * cp.p[end - 1] - cp.balance.debt, 0)

    # Check if cost of labor and investment can be financed from liquid assets
    NW_no_prod = (cp.balance.NW + cp.Dᵉ * cp.p[end] + cp.curracc.rev_dep 
                  - cp.debt_installments[1] - cp.balance.debt * globalparam.r)

    cp.possible_I = NW_no_prod + max_add_debt - TCLᵉ - TCE

    # If possible investments negative, decrease labor demand
    if cp.possible_I < 0

        cp.possible_I = 0

        poss_prod = (NW_no_prod + max_add_debt) / cop(cp.w̄[end], cp.π_LP, government.τᴱ, p_ep, cp.π_EE, government.τᶜ, cp.π_EF)
        poss_L = poss_prod / cp.π_LP
        cp.ΔLᵈ = poss_L - cp.L
    end
    
    # if NW_no_prod > TCLᵉ + TCE + cp.Iᵈ
    #     # All cost of costs can be paid from liquid assets. No additional debt needed.
    #     cp.cI = cp.Iᵈ

    # elseif NW_no_prod > TCLᵉ + TCE && NW_no_prod - TCLᵉ - TCE < cp.Iᵈ
    #     # Cost of labor can be paid from liquid assets, investment has to be partially 
    #     # funded from debt.
    #     cp.cI = NW_no_prod - TCLᵉ - TCE
    #     req_debt = cp.Iᵈ - (NW_no_prod - TCLᵉ - TCE)

    #     # Check if investment can be financed from additional debt, otherwise decrease investments
    #     if req_debt > max_add_debt
    #         if cp.EIᵈ > cp.cI + max_add_debt

    #             # Decrease amount of expansionary investment.
    #             poss_EI = cp.cI + max_add_debt
    #             cp.n_mach_ordered_EI = floor(Int64, cp.n_mach_ordered_EI * (poss_EI / cp.EIᵈ))
    #             cp.n_mach_ordered_RS = 0
    #             cp.mach_tb_repl = []

    #         else

    #             # Full expansion is possible, decrease amount of replacement investments
    #             poss_RS = cp.cI + max_add_debt - cp.EIᵈ
    #             cp.n_mach_ordered_RS = floor(Int64, cp.n_mach_ordered_RS * (poss_RS / cp.RSᵈ))

    #             if cp.n_mach_ordered_RS > 0
    #                 cp.mach_tb_repl = cp.mach_tb_repl[1:cp.n_mach_ordered_RS]
    #             else
    #                 cp.mach_tb_repl = []
    #             end

    #         end
    #     end

    # else
    #     # Cost of labor exceeds liquid assets. Check if enough additional debt available.
    #     # All investment cancelled.
    #     cp.cI = 0.0
    #     cp.n_mach_ordered_EI = 0
    #     cp.n_mach_ordered_RS = 0
    #     cp.mach_tb_repl = Vector{Machine}()

    #     if (NW_no_prod + max_add_debt) < (TCLᵉ + TCE)
    #         # Cost of labor exceeds expected liquid assets plus max additional debt. 
    #         # Decrease production quantity.
    #         # poss_prod = (NW_no_prod + max_add_debt) / (cp.w̄[end] / cp.π_LP + p_ep / cp.π_EE)
    #         poss_prod = (NW_no_prod + max_add_debt) / cop(cp.w̄[end], cp.π_LP, government.τᴱ, p_ep, cp.π_EE, government.τᶜ, cp.π_EF)
    #         poss_L = poss_prod / cp.π_LP
    #         cp.ΔLᵈ = poss_L - cp.L
    #     end
    # end

    # Based on final production decisions, update max offered wage
    update_wᴼ_max_cp!(cp)
end


"""
Lets cp make decision for kp out of available kp in brochures.
"""
function rank_producers_cp!(
    cp::ConsumerGoodProducer,
    government::Government, 
    b::Int64, 
    all_kp::Vector{Int64},
    ep::AbstractAgent,
    t::Int64,
    model::ABM
    )

    # In case of no brochures, pick a random kp
    if length(cp.kp_ids) == 0
        
        # Sample random kp id to add to known kp
        push!(cp.kp_ids, sample(all_kp))
    end

    all_cop = zeros(Float64, length(cp.kp_ids))
    for (i, kp_id) in enumerate(cp.kp_ids)

        brochure = get(model.kp_brochures, Symbol(kp_id), nothing)

        p_mach = get(brochure, :price, nothing)

        all_cop[i] = p_mach + b * cop(
                                cp.w̄[end], 
                                get(brochure, :A_LP, nothing), 
                                government.τᴱ, 
                                ep.p_ep[t], 
                                get(brochure, :A_EE, nothing), 
                                government.τᶜ, 
                                get(brochure, :A_EF, nothing)
                            )
    end


    # Sort kp ids based on cop
    cp.kp_ids .= sample(cp.kp_ids, Weights(1 ./ all_cop), length(cp.kp_ids); replace=false)
end

    # Choose kp based on brochures
    # all_cop .= (1 ./ all_cop .^ 2)
    # brochure = sample(cp.brochures, Weights(all_cop))
    # cp.kp_ranking[1] = brochure[:kp_id]
    # cp.kp_ranking = 

"""
Sort machines by cost of production
"""
function rank_machines_cp!(
    cp::ConsumerGoodProducer
)

    sort!(cp.Ξ , by = machine -> machine.cop; rev=true)
end


"""
Plans replacement investment based on age machines and available new machines
"""
function plan_replacement_cp!(
    cp::ConsumerGoodProducer,
    # kp_id::Int64,
    government::Government,
    globalparam::GlobalParam,
    ep::AbstractAgent,
    t::Int64,
    model::ABM
    )

    # desired_n_mach_replaced = 0
    # desired_n_mach_retired = 0
    n_machines_too_many = cp.n_machines > cp.Qᵉ ? cp.n_machines - cp.Qᵉ : 0

    # Get price and cost of production of chosen kp
    brochure = get(model.kp_brochures, Symbol(cp.kp_ids[1]), nothing)
    p_star = brochure[:price]
    c_star = cop(
                    cp.w̄[end], 
                    get(brochure, :A_LP, nothing), 
                    government.τᴱ, 
                    ep.p_ep[t], 
                    get(brochure, :A_EE, nothing), 
                    government.τᶜ, 
                    get(brochure, :A_EF, nothing)
                 )

    # Loop over machine stock, select which machines to replace
    for machine in cp.Ξ
        if machine.age >= globalparam.η
            # Machine has reached max age, decide if replaced or not
            if n_machines_too_many < machine.freq
                # No machines planned to be written off, replace old machine
                push!(cp.mach_tb_repl, machine)
                # desired_n_mach_replaced += 1
            else
                # Do not replace machine, let it be written off
                push!(cp.mach_tb_retired, machine)
                # desired_n_mach_retired == 1
                n_machines_too_many -= machine.freq
            end

        elseif (machine.cop != c_star && p_star / (machine.cop - c_star) <= globalparam.b 
                && machine.age > globalparam.b)
            # New machine cheaper to operate, replace old machine
            push!(cp.mach_tb_repl, machine)
            # desired_n_mach_replaced += 1
        end
    end

    cp.n_mach_ordered_RS = length(cp.mach_tb_repl)

    # return (min(max_mach_poss, desired_n_mach_replaced), desired_n_mach_retired)


    # If no new machines too expensive, no machines replaced
    # if cp.possible_I < p_star
    #     return (0, desired_n_mach_retired)
    # end

    # max_mach_poss = floor(Int64, cp.possible_I / p_star)

    # Sort to-be-replaced machines from lowest to highest production costs
    # sort!(cp.mach_tb_repl, by=machine->cop(cp.w̄[end], machine.A_LP, government.τᴱ, ep.p_ep[t], 
    #                                        machine.A_EE, government.τᶜ, machine.A_EF), rev=true)

    # Update total amount of to-be-replaces machines
    # add_mach_desired = length(cp.mach_tb_repl) - cp.n_mach_ordered_RS
    # println(cp.mach_tb_repl, " ", add_mach_desired, " ", max_mach_poss)
    # cp.n_mach_desired_RS = min(max_mach_poss, add_mach_desired)
    # c

    # @assert cp.n_mach_desired_RS >= 0

    # return (min(max_mach_poss, desired_n_mach_replaced), desired_n_mach_retired)

    # Compute investment amount corresponding to replacement investments
    # cp.RSᵈ = p_star * cp.n_mach_ordered_RS * globalparam.freq_per_machine
    # TODO replace with new computation
end


"""
Plans expansion investments based on expected production.
"""
function plan_expansion_cp!(
    cp::ConsumerGoodProducer,
    globalparam::GlobalParam,
    # roundnr::Int64,
    model::ABM
    )

    # If no more known capital producers, no orders possible
    # if roundnr > length(cp.kp_ids)
    #     cp.n_mach_desired_EI = 0
    #     return
    # end

    brochure = get(model.kp_brochures, Symbol(cp.kp_ids[1]), nothing)
    if cp.possible_I < brochure[:price]
        cp.n_mach_desired_EI = 0
        # cp.EIᵈ = 0.0
    end

    max_mach_poss = floor(Int64, cp.possible_I / brochure[:price])

    # if cp.Qᵉ > cp.n_machines && cp.cu > 0.8
    if cp.Qᵉ > cp.n_machines
        # cp.n_mach_ordered_EI = floor(Int64, (cp.Qᵉ - cp.n_machines) / globalparam.freq_per_machine)

        total_mach_desired = round(Int64, (cp.Qᵉ - cp.n_machines) / globalparam.freq_per_machine)
        add_mach_desired = total_mach_desired - cp.n_mach_ordered_EI
        cp.n_mach_desired_EI = max(min(max_mach_poss, add_mach_desired), 0)
        cp.n_mach_ordered_EI = cp.n_mach_desired_EI
        # TODO: CHANGE BACK THE DOUBLE VARS

        # println(cp.id, " ", roundnr, " ", add_mach_desired, " ", max_mach_poss, " ", total_mach_desired, " ", cp.n_mach_desired_RS, " ", cp.n_mach_ordered_EI, " ", cp.n_mach_ordered_RS)

        @assert cp.n_mach_desired_EI >= 0

        # cp.EIᵈ = brochure[:price] * cp.n_mach_ordered_EI * globalparam.freq_per_machine
        # TODO: Find a new way to compute this
    else
        cp.n_mach_desired_EI = 0
        # cp.EIᵈ = 0.0
    end
end


# """
# Plans replacement investment based on age machines and available new machines
# """
# function plan_replacement_cp!(
#     cp::ConsumerGoodProducer,
#     government::Government,
#     globalparam::GlobalParam,
#     ep::AbstractAgent,
#     # roundnr::Int64,
#     t::Int64,
#     model::ABM
#     )

#     # If no more known capital producers, no orders possible
#     if roundnr > length(cp.kp_ids)
#         cp.n_mach_desired_RS = 0
#         return
#     end

#     # Get price and cost of production of chosen kp
#     brochure = get(model.kp_brochures, Symbol(cp.kp_ids[roundnr]), nothing)
#     p_star = brochure[:price]
#     c_star = cop(
#                     cp.w̄[end], 
#                     get(brochure, :A_LP, nothing), 
#                     government.τᴱ, 
#                     ep.p_ep[t], 
#                     get(brochure, :A_EE, nothing), 
#                     government.τᶜ, 
#                     get(brochure, :A_EF, nothing)
#                  )

#     # See if machine stock too large in order to decide if need to be replaced
#     # n_machines_too_many = 0
#     # if cp.n_machines > cp.Qᵉ
#     #     n_machines_too_many = cp.n_machines - cp.Qᵉ
#     # end

#     n_machines_too_many = cp.n_machines > cp.Qᵉ ? cp.n_machines - cp.Qᵉ : 0

#     # Loop over machine stock, select which machines to replace
#     for machine in cp.Ξ

#         # If machine replaced by order in earlier round, do not check it
#         if machine ∈ cp.mach_tb_repl || machine ∈ cp.mach_tb_retired
#             continue
#         end

#         if machine.age >= globalparam.η
#             # Machine has reached max age, decide if replaced or not
#             if n_machines_too_many < machine.freq
#                 # No machines planned to be written off, replace old machine
#                 push!(cp.mach_tb_repl, machine)
#             else
#                 # Do not replace machine, let it be written off
#                 push!(cp.mach_tb_retired, machine)
#                 n_machines_too_many -= machine.freq
#             end

#         elseif (machine.cop != c_star && p_star / (machine.cop - c_star) <= globalparam.b 
#                 && machine.age > globalparam.b)
#             # New machine cheaper to operate, replace old machine
#             push!(cp.mach_tb_repl, machine)
#         end
#     end


#     # If no new machines too expensive, no machines replaced
#     if cp.possible_I < p_star
#         cp.n_mach_desired_RS = 0
#         return
#     end

#     max_mach_poss = floor(Int64, cp.possible_I / p_star)

#     # Sort to-be-replaced machines from lowest to highest production costs
#     sort!(cp.mach_tb_repl, by=machine->cop(cp.w̄[end], machine.A_LP, government.τᴱ, ep.p_ep[t], 
#                                            machine.A_EE, government.τᶜ, machine.A_EF), rev=true)

#     # Update total amount of to-be-replaces machines
#     add_mach_desired = length(cp.mach_tb_repl) - cp.n_mach_ordered_RS
#     # println(cp.mach_tb_repl, " ", add_mach_desired, " ", max_mach_poss)
#     cp.n_mach_desired_RS = min(max_mach_poss, add_mach_desired)

#     @assert cp.n_mach_desired_RS >= 0

#     # Compute investment amount corresponding to replacement investments
#     # cp.RSᵈ = p_star * cp.n_mach_ordered_RS * globalparam.freq_per_machine
#     # TODO replace with new computation
# end


# """
# Plans expansion investments based on expected production.
# """
# function plan_expansion_cp!(
#     cp::ConsumerGoodProducer,
#     globalparam::GlobalParam,
#     roundnr::Int64,
#     model::ABM
#     )

#     # If no more known capital producers, no orders possible
#     if roundnr > length(cp.kp_ids)
#         cp.n_mach_desired_EI = 0
#         return
#     end

#     brochure = get(model.kp_brochures, Symbol(cp.kp_ids[roundnr]), nothing)
#     if cp.possible_I < brochure[:price]
#         cp.n_mach_desired_EI = 0
#         # cp.EIᵈ = 0.0
#     end

#     max_mach_poss = floor(Int64, cp.possible_I / brochure[:price])

#     # if cp.Qᵉ > cp.n_machines && cp.cu > 0.8
#     if cp.Qᵉ > cp.n_machines
#         # cp.n_mach_ordered_EI = floor(Int64, (cp.Qᵉ - cp.n_machines) / globalparam.freq_per_machine)

#         total_mach_desired = round(Int64, (cp.Qᵉ - cp.n_machines) / globalparam.freq_per_machine)
#         add_mach_desired = total_mach_desired - cp.n_mach_ordered_EI
#         cp.n_mach_desired_EI = max(min(max_mach_poss, add_mach_desired), 0)
#         # println(cp.id, " ", roundnr, " ", add_mach_desired, " ", max_mach_poss, " ", total_mach_desired, " ", cp.n_mach_desired_RS, " ", cp.n_mach_ordered_EI, " ", cp.n_mach_ordered_RS)

#         @assert cp.n_mach_desired_EI >= 0

#         # cp.EIᵈ = brochure[:price] * cp.n_mach_ordered_EI * globalparam.freq_per_machine
#         # TODO: Find a new way to compute this
#     else
#         cp.n_mach_desired_EI = 0
#         # cp.EIᵈ = 0.0
#     end
# end


"""
Produces goods based on planned production and actual amount of hired workers
"""
function produce_goods_cp!(
    cp::ConsumerGoodProducer,
    ep::AbstractAgent,
    globalparam::GlobalParam,
    t::Int64
    )

    # If the cp does not need to use its complete capital stock, only use most productive 
    # machines
    n_machines_req = ceil(Int64, cp.Qˢ / globalparam.freq_per_machine)
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
    update_EU_TCE_cp!(cp, actual_π_EE, ep.p_ep[t])
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
    freq_per_machine::Int64,
    model::ABM
    )

    total_n_machines = cp.n_mach_ordered_EI + cp.n_mach_ordered_RS

    # Send orders for machines to kp
    if total_n_machines > 0 && hascapacity(model[cp.kp_ids[1]])
        receive_order_kp!(model[cp.kp_ids[1]], cp.id, total_n_machines, freq_per_machine)
    end
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
    t::Int64
    )

    cp.curracc.TCI += Iₜ

    # Replace old machines
    if length(cp.mach_tb_repl) > length(new_machines)
        # Not all to-be replaced machines were sent, only replace machines
        # that were delivered

        # Sort machines by cost of production, replace most expensive first
        sort!(cp.mach_tb_repl, by=machine->cp.w̄[end]/machine.A_LP + ep.p_ep[t]/machine.A_EE; rev=true)
        # filter!(machine -> machine ∉ cp.mach_tb_repl[1:length(new_machines)], cp.Ξ)
        for i in 1:length(new_machines)
            filter!(machine -> machine ≠ cp.mach_tb_repl[i], cp.Ξ)
        end
    else
        # All to-be replaced machines were sent, replace all machines and add
        # the additional machines as expansionary investment
        filter!(machine -> machine ∉ cp.mach_tb_repl, cp.Ξ)
    end

    append!(cp.Ξ, new_machines)
end


"""
Replaces cp, places cp in firm list of hh.
"""
function replace_bankrupt_cp!(
    bankrupt_cp::Vector{Int64},
    bankrupt_kp::Vector{Int64},
    all_hh::Vector{Int64},
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    globalparam::GlobalParam,
    indexfund::IndexFund,
    macro_struct::MacroEconomy,
    t::Int64,
    model::ABM
    )

    # Create vectors containing ids of non-bankrupt bp, lp and kp
    nonbankrupt_cp = setdiff(all_cp, bankrupt_cp)
    nonbankrupt_kp = setdiff(all_kp, bankrupt_kp)

    # Compute average number of machines and NW for non-bankrupt cp
    avg_n_machines = mean(cp_id -> model[cp_id].n_machines, nonbankrupt_cp)
    avg_NW = mean(cp_id -> model[cp_id].balance.NW, nonbankrupt_cp)

    # Make weights for allocating cp to hh
    # Minimum is taken to avoid weird outcomes when all bp and lp went bankrupt
    weights_hh_cp = map(hh_id -> min(1, 1 / length(model[hh_id].cp)), all_hh)
    weights_kp = map(kp_id -> max(model[kp_id].f[end], 0.01), nonbankrupt_kp)

    n_bankrupt_cp = length(bankrupt_cp)

    # Sample all NW coefficients and capital coefficients
    capital_coefficients = rand(Uniform(globalparam.φ1, globalparam.φ2), n_bankrupt_cp)
    NW_coefficients = rand(Uniform(globalparam.φ3, globalparam.φ4), n_bankrupt_cp)

    # New cp receive an advanced type of machine, first select kp ids proportional
    # to their market share. cp can also select kp ids that went bankrupt in this 
    # period, as these producers have already been replaced with new companies
    kp_choice_ids = zeros(Int64, n_bankrupt_cp)
    kp_choice_ps = zeros(Float64, n_bankrupt_cp)
    all_n_machines = zeros(Int64, n_bankrupt_cp)

    for i in 1:n_bankrupt_cp
        # Decide from which kp to buy
        kp_choice_ids[i] = sample(all_kp, Weights(weights_kp))
        kp_choice_ps[i] = model[kp_choice_ids[i]].p[end]

        # Compute the number of machines each cp will buy
        all_n_machines[i] = floor(Int64, capital_coefficients[i] * avg_n_machines / globalparam.freq_per_machine)
    end

    # Compute share of investments that can be paid from the investment fund
    req_NW = (avg_NW .* NW_coefficients) .+ (all_n_machines .* (kp_choice_ps .* globalparam.freq_per_machine))
    all_req_NW = sum(req_NW)
    frac_NW_if = decide_investments_if!(indexfund, all_req_NW, t)

    n_kp_sample = min(length(weights_kp), 10)

    for (cp_i, cp_id) in enumerate(bankrupt_cp)

        # Sample what the size of the capital stock will be
        D = macro_struct.cu[t] * all_n_machines[cp_i] * globalparam.freq_per_machine

        # In the first period, the cp has no machines yet, these are delivered at the end
        # of the first period
        new_cp = initialize_cp(
                    cp_id,
                    cp_i,
                    # t + 1,
                    Vector{Machine}(),
                    model;
                    D=D,
                    w=macro_struct.w_avg[t],
                    f=0.0
                )

        # Order machines at kp of choice
        new_cp.kp_ids = sample(nonbankrupt_kp, Weights(weights_kp), n_kp_sample; replace=false)
        new_cp.n_mach_desired_EI = all_n_machines[cp_i]

        update_wᴼ_max_cp!(new_cp)

        # Augment the balance with acquired NW and K
        new_cp.balance.NW = req_NW[cp_i]

        # Borrow remaining required funds for the machine, the other part of the 
        # funds come from the investment fund
        borrow_funds_p!(new_cp, (1 - frac_NW_if) * req_NW[cp_i], globalparam.b)

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

    cp.Dᵉ = cp.age > 1 ? ω * cp.Dᵉ + (1 - ω) * (cp.D[end] + cp.Dᵁ[end]) : cp.Dᵉ
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

    if length(cp.Ξ) > 0
        cp.Qᵉ = ω * cp.Qᵉ + (1 - ω) * ((1 + ι) * cp.Dᵉ)
    end
end


"""
Updates capital stock n_machines
"""
function update_n_machines_cp!(
    cp::ConsumerGoodProducer,
    freq_per_machine::Int64
    )

    # Retire old machines that are not replaced
    if length(cp.mach_tb_retired) > 0
        setdiff!(cp.Ξ, cp.mach_tb_retired)
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


# """
# Computes the markup rate μ based on the market share f.
# """
# function update_μ_cp!(
#     cp::ConsumerGoodProducer,
#     ϵ_μ::Float64,
#     t::Int64
#     )

#     new_μ = cp.μ[end]
#     shock = ϵ_μ * rand()

#     if cp.age > 2 && t > 2
#         dp = cp.μ[end] - cp.μ[end - 1]
#         dΠ = cp.Π[end] - cp.Π[end - 1]
#         new_μ *= (1 + shock * sign(dp) * sign(dΠ))
#     else
#         new_μ *= (1 + shock * sample([-1., 1.]))
#     end

#     shift_and_append!(cp.μ, new_μ)
# end


"""
Compute production cost per unit c
"""
function compute_c_cp!(
    cp::ConsumerGoodProducer,
    p_ep::Float64,
    τᴱ::Float64,
    τᶜ::Float64
    )

    if cp.L + cp.ΔLᵈ > 0
        c = cop(cp.w̄[end], cp.π_LP, τᴱ, p_ep, cp.π_EE, τᶜ, cp.π_EF)
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
Computes the desired labor supply Lᵈ and the change in labor supply 
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
    
    cp.wᴼ_max = cp.π_LP * cp.p[end]
end


"""
Updates energy use for production
"""
function update_EU_TCE_cp!(
    cp::ConsumerGoodProducer,
    actual_π_EE::Float64,
    p_ep::Float64
    )

    cp.EU = length(cp.Ξ) > 0 ? cp.Q[end] / actual_π_EE : 0.0
    cp.curracc.TCE = p_ep * cp.EU
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


function checkupdateprice(
    cp_id::Int64, 
    t::Int64,
    n_hh::Int64, 
    n_cp::Int64,
    p_time_rigid::Int64
)

    id_rescaled = (cp_id - n_hh) ÷ (n_cp / p_time_rigid) + 1
    t_rescaled = t % p_time_rigid

    return id_rescaled == t_rescaled
end


function add_kp_cp!(
    cp::ConsumerGoodProducer, 
    kp_id::Int64
)

    if kp_id ∉ cp.kp_ids
        push!(cp.kp_ids, kp_id)
    end
end


"""
Sets number of ordered and desired machines to zero
"""
function reset_desired_ordered_machines_cp!(
    cp::ConsumerGoodProducer
)

    cp.n_mach_ordered_EI = 0
    cp.n_mach_desired_EI = 0
    cp.n_mach_ordered_RS = 0
    cp.n_mach_desired_RS = 0

    cp.mach_tb_repl = []
    cp.mach_tb_retired = []
end


"""

"""
function process_machine_order_cp!(
    cp::ConsumerGoodProducer, 
    order::Int64,
    price::Float64
)

    # First satisfy replacements
    cp.n_mach_ordered_RS += min(cp.n_mach_desired_RS, order)
    cp.mach_tb_repl = cp.mach_tb_repl[1:cp.n_mach_ordered_RS]

    # Secondly satisfy expansions
    cp.n_mach_ordered_EI += max(order - cp.n_mach_desired_RS, 0)

    # println(cp.n_mach_ordered_EI, " ", order, " ", cp.n_mach_ordered_RS)

    @assert cp.n_mach_ordered_EI >= 0

    # Decrease possible investments
    cp.possible_I -= price * order

    # Set desired number of machines back to zero for next round
    cp.n_mach_desired_RS = 0
    cp.n_mach_desired_EI = 0
end


function update_cop_machines_cp!(
    cp::ConsumerGoodProducer, 
    government::Government, 
    ep::AbstractAgent,
    t::Int64
)

    for machine in cp.Ξ
        update_cop_machine!(machine, cp, government, ep, t)
    end
end