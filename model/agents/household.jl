# @pyimport scipy.stats as stats

mutable struct Household <: AbstractAgent

    id :: Int                   # global id

    # Employment variables
    employed :: Bool            # is employed
    employer_id :: Int          # id of employer
    L :: Float64                # labor units in household
    w :: Vector{Float64}       # wage
    wˢ :: Float64               # satisfying wage
    wʳ :: Float64               # requested wage
    T_unemp :: Int              # time periods unemployed

    # Income and wealth variables
    I :: Vector{Float64}        # hist income
    Iᵀ :: Vector{Float64}       # hist taxed income
    # S :: Vector{Float64}        # total savings
    s :: Float64                # savings rate
    W :: Vector{Float64}        # wealth or cash on hand
    # Wʳ :: Vector{Float64}       # real wealth or cash on hand

    # Consumption variables
    C :: Vector{Float64}        # budget
    # N_goods :: Float64          # number of goods bought
    bp :: Vector{Int}           # connected cp basic goods
    lp :: Vector{Int}           # connected cp luxury goods
    unsat_dem :: Vector         # unsatisfied demands
    P̄ :: Float64                # weighted average price of bp
    c_L :: Float64              # share of income used to buy luxury goods

end

function initialize_hh(
    id::Int,
    )::Household

    hh = Household(
        id,                     # global id

        false,                  # bool: employed
        0,                      # id of employer
        100,                    # L: labor units
        ones(4),                # w: wage
        1.0,                    # wˢ: satisfying wage
        1.0,                    # wʳ: requested wage
        0,                      # T_unemp: time periods unemployed

        [100],                  # I: hist income
        [],                     # Iᵀ: hist taxed income
        # [10],                   # S: total savings
        0,                      # s: savings rate
        [50],                  # W: wealth or cash on hand
        # [100],                  # Wʳ: real wealth or cash on hand

        [],                     # C: budget
        # 0.0,                    # N_goods: number of goods bought
        Vector{Int}(),          # all_cp_B: connected cp basic goods
        Vector{Int}(),          # all_cp_L: connected cp luxury goods
        Vector(),
        0,                      # P̄: weighted average price of bp
        0.5,                    # c_L: share of budget used to buy luxury goods
    )
    return hh
end


"""
Uniformly samples basic and luxury good producers to be in trading network.
"""
function select_bp_lp_hh!(
    hh::Household,
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    n_bp::Int,
    n_lp::Int
    )

    selected_bp = sample(all_bp, n_bp)
    hh.bp = selected_bp

    selected_lp = sample(all_lp, n_lp)
    hh.lp = selected_lp
end


# """
# Updates current wealth level W, equal to the current cash on hand
# """
# function update_W_hh!(
#     hh::Household
#     )
    
#     W = hh.Iᵀ[end] + hh.S[end]
#     push!(hh.W, W)
# end


"""
Sets consumption budget based on current wealth level
"""
function set_consumption_budget_hh!(
    hh::Household,
    all_W_hh::Vector{Float64},
    global_param::GlobalParam,
    model::ABM
    )

    # Update average price level of bp and lp
    update_average_price_hh!(hh, model)

    # Update real wealth level
    # update_real_wealth_hh!(hh)
    # TODO: check if this is still needed

    # Update share of goods to bg and lg
    update_share_goodtypes_hh!(hh, global_param.c_L_max, global_param.a_σ, global_param.b_σ)

    # Compute consumption budget
    compute_consumption_budget_hh!(hh, global_param.α_cp, all_W_hh)
end



"""
Computes average price level of available bp and lp
"""
function update_average_price_hh!(
    hh::Household,
    model::ABM
    )

    P̄_bp = mean(bp_id -> model[bp_id].p[end], hh.bp)
    P̄_lp = mean(lp_id -> model[lp_id].p[end], hh.lp)
    hh.P̄ = hh.c_L * P̄_lp + (1 - hh.c_L) * P̄_bp
end


# """
# Updates real wealth level based on price division of last period
# """
# function update_real_wealth_hh!(
#     hh::Household
#     )

#     Wʳ = hh.W[end] / hh.P̄[end]
#     push!(hh.Wʳ, Wʳ)
# end


"""
Defines logistic function that determines share of each good type
"""
function update_share_goodtypes_hh!(
    hh::Household,
    c_L_max::Float64,
    a_σ::Float64,
    b_σ::Float64
    )

    # TODO: change this back once switching is implemented
    # hh.c_L = c_L_max / (1 + exp(-(hh.Wʳ[end]/a_σ - b_σ)))
    hh.c_L = 0.5
end


"""
Computes consumption budget, updates savings rate
"""
function compute_consumption_budget_hh!(
    hh::Household,
    α_cp::Float64,
    all_W_hh::Vector{Float64}
    )

    if hh.W[end] > 0 && hh.Iᵀ[end] > 0
        # percentile = stats.percentileofscore(all_W_hh, hh.W[end])
        # percentile = 50
        # frac_cons = percentile^α_cp / percentile
        # println("$frac_cons, $percentile")
        C = min(hh.P̄[end] * (hh.W[end] / hh.P̄[end])^α_cp, hh.W[end])
        # C = hh.Iᵀ[end]
        # C = frac_cons * hh.W[end] / 100
        s = (hh.Iᵀ[end] - C) / hh.Iᵀ[end]
    else
        C = 0
        s = 0
    end

    push!(hh.C, C)
    hh.s = s
end


"""
Places orders at bp and lp
"""
function place_orders_hh!(
    hh_p::Vector{Int},
    C_i::Float64,
    cp_inventories,
    p_with_inventory::Vector{Int},
    global_param::GlobalParam,
    model::ABM
    )

    # If none of the known suppliers has products in stock, randomly select
    # other suppliers until demand can be met.
    poss_p = Dict(p_id => 1 / model[p_id].p[end]^2 for p_id ∈ intersect(Set(hh_p), Set(p_with_inventory)))

    add_p_id = 0

    # As long as the current producers do not have enough inventory and there are still
    # possible producers to sample, randomly sample producers and add to pool of possible cp
    while (sum(map(p_id -> cp_inventories[p_id], collect(keys(poss_p)))) < C_i && 
           length(poss_p) != length(p_with_inventory))

        add_p_id = sample(setdiff(p_with_inventory, keys(poss_p)))
        poss_p[add_p_id] = 1 / model[add_p_id].p[end] 
    end

    # Place orders based on the availability of goods
    chosen_amount = 0
    chosen_p_id = 0
    C_per_day = C_i / global_param.n_cons_market_days

    p_orders = Dict()
    
    while C_i > 0 && length(poss_p) > 0
        chosen_p_id = sample(collect(keys(poss_p)), Weights(collect(values(poss_p))))
        chosen_amount = min(min(C_i, C_per_day), cp_inventories[chosen_p_id])

        if chosen_p_id ∉ keys(p_orders)
            p_orders[chosen_p_id] = chosen_amount
        else
            p_orders[chosen_p_id] += chosen_amount
        end

        C_i -= chosen_amount
        cp_inventories[chosen_p_id] -= chosen_amount

        if cp_inventories[chosen_p_id] == 0
            delete!(poss_p, chosen_p_id)
            filter!(p_id -> p_id ≠ chosen_p_id, p_with_inventory)
        end
    end

    return p_orders, cp_inventories, p_with_inventory
end


"""
Household receives ordered cg and mutates balance
"""
function receive_order_hh!(
    hh::Household,
    cp_id::Int,
    tot_price::Float64,
    share_fulfilled::Float64
    )

    # Decrease wealth with total price paid
    hh.W[end] -= tot_price

    # If full demand not fulfilled, add cp to unsatisfied demand
    if share_fulfilled < 1.0
        push!(hh.unsat_dem, (cp_id, 1 - share_fulfilled))
    end
end


"""
Updates satisfying wage wˢ and requested wage wʳ
"""
function update_sat_req_wage_hh!(
    hh::Household, 
    ϵ::Float64, 
    UB::Float64
    )

    # T = 4
    # if length(hh.w) > T
    #     hh.wˢ = mean(hh.w[end-T:end])/hh.L
    # else
    #     hh.wˢ = mean(hh.w)/hh.L
    # end

    # Try to use adaptive wˢ
    ωwˢ = 0.8
    hh.wˢ = ωwˢ * hh.wˢ + (1 - ωwˢ) * hh.I[end] / hh.L
    # hh.wˢ = ωwˢ * hh.wˢ + (1 - ωwˢ) * hh.w[end]

    if hh.employed
        hh.wʳ = hh.w[end] * (1 + ϵ)
    else
        hh.wʳ = max(UB/hh.L, hh.wˢ)
    end
end


"""
Lets households get income, either from UB or wage
"""
function get_income_hh!(
    hh::Household, 
    amount::Float64
    )
    # if isnan(amount)
    #     println("yeet ", hh.employer_id)
    # end
    # println(amount)
    push!(hh.I, amount)
    if hh.employed
        hh.w[1:3] = hh.w[2:4]
    end
    # println("1 ", amount, " ", hh.W[end])
    # push!(hh.W, amount + hh.W[end])
    # println("2 ", amount, " ", hh.W[end])
end


"""
Updates household wealth
"""
function update_wealth_hh!(
    hh::Household
    )

    W = hh.W[end] + hh.Iᵀ[end]
    push!(hh.W, W)
end


"""
Sets household to be unemployed.
"""
function set_unemployed_hh!(
    hh::Household
    )

    hh.employed = false
    hh.employer_id = 0
end


"""
Lets employee be hired when previously unemployed, saves employer id and new earned wage.
"""
function set_employed_hh!(
    hh::Household, 
    wᴼ::Float64,
    employer_id::Int,
    )

    hh.employed = true
    hh.employer_id = employer_id
    hh.T_unemp = 0
    # push!(hh.w, wᴼ)
    hh.w[1:3] = hh.w[2:4]
    hh.w[4] = wᴼ
end


"""
Changes employer for households that were already employed.
"""
function change_employer_hh!(
    hh::Household,
    wᴼ::Float64,
    employer_id::Int
    )

    hh.employer_id = employer_id
    # push!(hh.w, wᴼ)
    hh.w[1:3] = hh.w[2:4]
    hh.w[4] = wᴼ
end


"""
Removes bankrupt producers from set of producers.
"""
function remove_bankrupt_producers_hh!(
    hh::Household,
    bankrupt_bp::Vector{Int},
    bankrupt_lp::Vector{Int}
    )

    filter!(bp_id -> bp_id ∉ bankrupt_bp, hh.bp)
    filter!(lp_id -> lp_id ∉ bankrupt_lp, hh.lp)
end


"""
Decides whether to switch to other cp
"""
function decide_switching_hh!(
    hh::Household,
    ψ_Q::Float64,
    ψ_P::Float64,
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    model::ABM
    )

    # Check if demand was constrained and for chance of changing cp
    if length(hh.unsat_dem) > 0 && rand() < ψ_Q

        # Pick a supplier to change, first set up weights inversely proportional
        # to supplied share of goods
        weights = map(p -> 1/max(p[2], 0.001), hh.unsat_dem)

        # Sample producer to replace
        p_id_replaced = sample(hh.unsat_dem, Weights(weights))[1]

        # Check if replaced supplier is bp or lp, sample new supplier in correct
        # category and replace in set of hh suppliers.
        if p_id_replaced in hh.bp
            filter!(p_id -> p_id ≠ p_id_replaced, hh.bp)

            # Add new bp if list not already too long
            if length(hh.bp) < 10
                p_id_new = sample(setdiff(all_bp, hh.bp))
                push!(hh.bp, p_id_new)
            end
        else
            filter!(p_id -> p_id ≠ p_id_replaced, hh.lp)

            # Add new bp if list not already too long
            if length(hh.lp) < 10
                p_id_new = sample(setdiff(all_lp, hh.lp))
                push!(hh.lp, p_id_new)
            end
        end
    end

    n_attempts = 5

    # Check if household will look for a better price
    if rand() < ψ_P

        # Randomly select a supplier that may be replaced
        # p_id_candidate1 = max(p_id -> model[p_id].p[end], vcat(hh.bp, hh.lp))
        p_id_candidate1 = sort(vcat(hh.bp, hh.lp), by = p_id -> model[p_id].p[end])[end]

        # Randomly pick another candidate from same type and see if price is lower
        if p_id_candidate1 ∈ hh.bp
            # TODO make this weighted (if needed)
            # TODO see if you dont always want to sample from already known producers
            p_id_candidate2 = sort(sample(setdiff(all_bp, hh.bp), n_attempts), 
                                   by = bp_id -> model[bp_id].p[end])[1]
            
            # Replace supplier if price of other supplier is lower 
            if model[p_id_candidate2].p[end] < model[p_id_candidate1].p[end]
                hh.bp[findall(x->x==p_id_candidate1, hh.bp)] .= p_id_candidate2
                # filter!(p_id -> p_id ≠ p_id_candidate1, hh.bp)
                # push!(hh.bp, p_id_candidate2)
            end
        else
            p_id_candidate2 = sort(sample(setdiff(all_lp, hh.lp), n_attempts), 
                                   by = lp_id -> model[lp_id].p[end])[1]
        
            # Replace supplier if price of other supplier is lower
            if model[p_id_candidate2].p[end] < model[p_id_candidate1].p[end]
                hh.lp[findall(x->x==p_id_candidate1, hh.lp)] .= p_id_candidate2
                # filter!(p_id -> p_id ≠ p_id_candidate1, hh.lp)
                # push!(hh.lp, p_id_candidate2)
            end
        end
    end
end


"""
Refills amount of bp and lp in amount is below minimum. Randomly draws suppliers
    inversely proportional to prices.
"""
function refill_suppliers_all_hh!(
    all_hh::Vector{Int},
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    model::ABM;
    n_bp_hh=7::Int,
    n_lp_hh=7::Int
    )

    # Check amount of bp and lp of all hh, if insufficient replenish by randomly
    # drawing other suppliers based on weights
    for hh_id in all_hh

        # Check bp, sample if insufficient
        if length(model[hh_id].bp) < n_bp_hh

            # Determine which bp are available
            n_add_bp = n_bp_hh - length(model[hh_id].bp)
            poss_bp = filter(p_id -> p_id ∉ model[hh_id].bp, all_bp)

            # Determine weights based on prices, sample and add
            weights_bp = map(bp_id -> 1 / model[bp_id].p[end], poss_bp)
            add_bp = sample(poss_bp, Weights(weights_bp), n_add_bp)
            model[hh_id].bp = vcat(model[hh_id].bp, add_bp)
        end

        # Check lp, sample if insufficient
        if length(model[hh_id].lp) < n_lp_hh

            # Determine which lp are available
            n_add_lp = n_lp_hh - length(model[hh_id].lp)
            poss_lp = filter(p_id -> p_id ∉ model[hh_id].lp, all_lp)

            # Determine weights based on prices, sample and add
            weights_lp = map(lp_id -> 1 / model[lp_id].p[end], poss_lp)
            add_lp = sample(poss_lp, Weights(weights_lp), n_add_lp)
            model[hh_id].lp = vcat(model[hh_id].lp, add_lp)
        end
    end
end