@with_kw mutable struct Household <: AbstractAgent

    id :: Int                                   # global id

    # Employment variables
    employed::Bool = false                      # is employed
    employer_id::Union{Int} = 0                 # id of employer
    L::Float64 = 100.0                          # labor units in household
    w::Vector{Float64} = ones(Float64, 4)       # wage
    wˢ::Float64 = 1.0                           # satisfying wage
    wʳ::Float64 = 1.0                           # requested wage
    T_unemp::Int = 0                            # time periods unemployed

    # Income and wealth variables
    I::Float64 = 100                            # hist income
    Iᵀ::Float64 = 100    # change according to tax rate      # hist taxed income
    s::Float64 = 0.0                            # savings rate
    W::Float64 = 50                             # wealth or cash on hand

    # Consumption variables
    C::Float64 = 0.0                           # budget
    bp::Vector{Int} = Vector{Int}()            # connected cp basic goods
    lp::Vector{Int} = Vector{Int}()            # connected cp luxury goods
    unsat_dem::Vector = Vector()               # unsatisfied demands
    P̄::Float64 = 1.0                           # weighted average price of bp
    c_L::Float64 = 0.5                         # share of income used to buy luxury goods
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

    hh.bp = sample(all_bp, n_bp)
    hh.lp = sample(all_lp, n_lp)
end


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

    if hh.W > 0
        hh.C = min(hh.P̄[end] * (hh.W / hh.P̄[end])^α_cp, hh.W)
        # hh.C = hh.W
        hh.s = hh.Iᵀ > 0 ? (hh.Iᵀ - hh.C) / hh.Iᵀ : -1.0
    else
        hh.C = 0.0
        hh.s = 0.0
    end
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
    model::ABM,
    to
    )

    # If none of the known suppliers has products in stock, randomly select
    # other suppliers until demand can be met.

    # Ugly selection to boost performance.
    poss_p = Dict(p_id => 1 / model[p_id].p[end]^2 for p_id ∈ hh_p)
    for p_id in hh_p
        if p_id ∉ p_with_inventory
            delete!(poss_p, p_id)
        end
    end

    add_p_id = 0
    poss_p_ids = collect(keys(poss_p))
    sum_poss_p = length(poss_p_ids) > 0 ? sum(p_id -> cp_inventories[p_id], poss_p_ids) : 0

    # As long as the current producers do not have enough inventory and there are still
    # possible producers to sample, randomly sample producers and add to pool of possible cp
    while sum_poss_p < C_i && length(poss_p) != length(p_with_inventory)

        add_p_id = sample(p_with_inventory)
        while add_p_id ∈ poss_p_ids
            add_p_id = sample(p_with_inventory)
        end 

        poss_p[add_p_id] = 1 / model[add_p_id].p[end]
        sum_poss_p += cp_inventories[add_p_id]
        push!(poss_p_ids, add_p_id) 
    end

    # Place orders based on the availability of goods
    chosen_amount = 0
    chosen_p_id = 0
    C_per_day = C_i / global_param.n_cons_market_days

    p_orders = Dict(cp_id => 0.0 for cp_id in poss_p_ids)
    weights = collect(values(poss_p))
    
    n_round = 1
    n_round_stop::Int = global_param.n_cons_market_days * 1.5

    while C_i > 0 && length(poss_p) > 0 && n_round < n_round_stop
        chosen_p_id = sample(poss_p_ids, Weights(weights))
        chosen_amount = min(C_i, C_per_day, cp_inventories[chosen_p_id])
        p_orders[chosen_p_id] += chosen_amount

        C_i -= chosen_amount
        cp_inventories[chosen_p_id] -= chosen_amount

        if cp_inventories[chosen_p_id] == 0
            filter!(p_id -> p_id ≠ chosen_p_id, poss_p_ids)
            filter!(weight -> weight ≠ poss_p[chosen_p_id], weights)
            delete!(poss_p, chosen_p_id)
            filter!(p_id -> p_id ≠ chosen_p_id, p_with_inventory)
        end
        n_round += 1
    end

    return p_orders, cp_inventories, p_with_inventory
end


"""
Household receives ordered cg and mutates balance
"""
function receive_ordered_goods_hh!(
    hh::Household,
    cp_id::Int,
    tot_price::Float64,
    share_fulfilled::Float64
    )

    # Decrease wealth with total price paid
    hh.W -= tot_price

    # If full demand not fulfilled, add cp to unsatisfied demand
    if ceil(share_fulfilled; digits=3) < 1.0
        println(share_fulfilled)
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
    hh.wˢ = ωwˢ * hh.wˢ + (1 - ωwˢ) * hh.Iᵀ / hh.L

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

    hh.I = amount
    if hh.employed
        shift_and_append!(hh.w, hh.w[end])
    end
end


"""
Updates household wealth
"""
function update_wealth_hh!(
    hh::Household
    )

    hh.W += hh.Iᵀ
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
    shift_and_append!(hh.w, wᴼ)
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
    shift_and_append!(hh.w, wᴼ)
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
function decide_switching_all_hh!(
    global_param::GlobalParam,
    all_hh::Vector{Int},
    all_p::Vector{Int},
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    model::ABM
    )

    # all_P_weights = Dict(p_id => 1/max(model[p_id].p[end], 0.001) for p_id in all_p)

    for hh_id in all_hh
        # Check if demand was constrained and for chance of changing cp
        if length(model[hh_id].unsat_dem) > 0 && rand() < global_param.ψ_Q

            # Pick a supplier to change, first set up weights inversely proportional
            # to supplied share of goods
            weights = map(p -> 1/max(p[2], 0.001), model[hh_id].unsat_dem)

            # Sample producer to replace
            p_id_replaced = sample(model[hh_id].unsat_dem, Weights(weights))[1]

            # Check if replaced supplier is bp or lp, sample new supplier in correct
            # category and replace in set of hh suppliers.
            if p_id_replaced ∈ model[hh_id].bp
                filter!(p_id -> p_id ≠ p_id_replaced, model[hh_id].bp)

                # Add new bp if list not already too long
                if length(model[hh_id].bp) < 10
                    # p_id_new = sample(setdiff(all_bp, model[hh_id].bp))
                    p_id_new = sample(all_bp)
                    push!(model[hh_id].bp, p_id_new)
                end
            else
                filter!(p_id -> p_id ≠ p_id_replaced, model[hh_id].lp)

                # Add new bp if list not already too long
                if length(model[hh_id].lp) < 10
                    # p_id_new = sample(setdiff(all_lp, model[hh_id].lp))
                    p_id_new = sample(all_lp)
                    push!(model[hh_id].lp, p_id_new)
                end
            end
        end

        # Check if household will look for a better price
        if rand() < global_param.ψ_P

            # Randomly select a supplier that may be replaced
            # p_id_candidate1 = sort(vcat(model[hh_id].bp, model[hh_id].lp), by=p_id -> model[p_id].p[end])[end]
            p_id_candidate1 = sample(vcat(model[hh_id].bp, model[hh_id].lp))

            # Randomly pick another candidate from same type and see if price is lower
            if p_id_candidate1 ∈ model[hh_id].bp
                # TODO make this weighted (if needed)
                # TODO see if you dont always want to sample from already known producers

                # Ugly sample to boost performance
                p_id_candidate2 = sample(all_bp)
                while p_id_candidate2 ∈ model[hh_id].bp
                    p_id_candidate2 = sample(all_bp)
                end
                
                # Replace supplier if price of other supplier is lower 
                if model[p_id_candidate2].p[end] < model[p_id_candidate1].p[end]
                    model[hh_id].bp[findfirst(x->x==p_id_candidate1, model[hh_id].bp)] = p_id_candidate2
                end
            else
                p_id_candidate2 = sample(all_lp)
                while p_id_candidate2 ∈ model[hh_id].lp
                    p_id_candidate2 = sample(all_lp)
                end
            
                # Replace supplier if price of other supplier is lower
                if model[p_id_candidate2].p[end] < model[p_id_candidate1].p[end]
                    model[hh_id].lp[findfirst(x->x==p_id_candidate1, model[hh_id].lp)] = p_id_candidate2
                end
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