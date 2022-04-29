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
    cp::Vector{Int} = Vector{Int}()            # connected cp
    unsat_dem::Vector = Vector()               # unsatisfied demands
    P̄::Float64 = 1.0                           # weighted average price of bp
    c_L::Float64 = 0.5                         # share of income used to buy luxury goods
end


"""
Uniformly samples cp to be in trading network.
"""
function select_cp_hh!(
    hh::Household,
    all_cp::Vector{Int},
    n_cp_hh::Int
    )

    hh.cp = sample(all_cp, n_cp_hh)
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

    hh.P̄ = mean(cp_id -> model[cp_id].p[end], hh.cp)
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
        hh.s = hh.Iᵀ > 0 ? (hh.Iᵀ - hh.C) / hh.Iᵀ : -1.0
    else
        hh.C = 0.0
        hh.s = 0.0
    end
end


"""
Places orders at cp
"""
function place_orders_hh!(
    hh_p::Vector{Int},
    C::Float64,
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
    while sum_poss_p < C && length(poss_p) != length(p_with_inventory)

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
    C_per_day = C / global_param.n_cons_market_days

    p_orders = Dict(cp_id => 0.0 for cp_id in poss_p_ids)
    weights = collect(values(poss_p))
    
    n_round = 1
    n_round_stop::Int = global_param.n_cons_market_days * 1.5

    while C > 0 && length(poss_p) > 0 && n_round < n_round_stop
        chosen_p_id = sample(poss_p_ids, Weights(weights))
        chosen_amount = min(C, C_per_day, cp_inventories[chosen_p_id])
        p_orders[chosen_p_id] += chosen_amount

        C -= chosen_amount
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


# """
# Household receives ordered cg and mutates balance
# """
# function receive_ordered_goods_hh!(
#     hh::Household,
#     cp_id::Int,
#     tot_price::Float64,
#     share_fulfilled::Float64
#     )

#     # Decrease wealth with total price paid
#     hh.W -= tot_price

#     # If full demand not fulfilled, add cp to unsatisfied demand
#     if ceil(share_fulfilled; digits=3) < 1.0
#         println(share_fulfilled)
#         push!(hh.unsat_dem, (cp_id, 1 - share_fulfilled))
#     end
# end

"""
Household receives ordered cg and mutates balance
"""
function receive_ordered_goods_hh!(
    hh::Household,
    tot_sales::Float64,
    unsat_demand::Vector{Float64},
    hh_D::Vector{Float64},
    all_cp::Vector{Int}
    # cp_id::Int,
    # tot_price::Float64,
    # share_fulfilled::Float64
    )

    # Decrease wealth with total sold goods
    hh.W -= tot_sales

    for i in findall(cp_id -> cp_id ∈ hh.cp, all_cp)
        if unsat_demand[i] > 0.0
            push!(hh.unsat_dem, (all_cp[i], unsat_demand[i] / hh_D[i]))
        end
    end

    # # If full demand not fulfilled, add cp to unsatisfied demand
    # if ceil(share_fulfilled; digits=3) < 1.0
    #     println(share_fulfilled)
    #     push!(hh.unsat_dem, (cp_id, 1 - share_fulfilled))
    # end
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
    bankrupt_cp::Vector{Int}
    )

    filter!(cp_id -> cp_id ∉ bankrupt_cp, hh.cp)
end


"""
Decides whether to switch to other cp
"""
function decide_switching_all_hh!(
    global_param::GlobalParam,
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    all_p::Vector{Int},
    n_cp_hh::Int,
    model::ABM
    )

    for hh_id in all_hh
        # Check if demand was constrained and for chance of changing cp
        if length(model[hh_id].unsat_dem) > 0 && rand() < global_param.ψ_Q

            # Pick a supplier to change, first set up weights inversely proportional
            # to supplied share of goods
            weights = map(p -> 1/max(p[2], 0.001), model[hh_id].unsat_dem)

            # Sample producer to replace
            p_id_replaced = sample(model[hh_id].unsat_dem, Weights(weights))[1]

            filter!(p_id -> p_id ≠ p_id_replaced, model[hh_id].cp)

            # Add new cp if list not already too long
            if length(model[hh_id].cp) < n_cp_hh
                
                p_id_new = sample(all_cp)
                while p_id_new ∈ model[hh_id].cp
                    p_id_new = sample(all_cp)
                end
                push!(model[hh_id].cp, p_id_new)
            end

        end

        # Check if household will look for a better price
        if rand() < global_param.ψ_P

            # Randomly select a supplier that may be replaced
            p_id_candidate1 = sample(model[hh_id].cp)

            # Randomly pick another candidate from same type and see if price is lower
            # Ugly sample to boost performance
            p_id_candidate2 = sample(all_cp)
            while p_id_candidate2 ∈ model[hh_id].cp
                p_id_candidate2 = sample(all_cp)
            end
            
            # Replace supplier if price of other supplier is lower 
            if model[p_id_candidate2].p[end] < model[p_id_candidate1].p[end]
                model[hh_id].cp[findfirst(x->x==p_id_candidate1, model[hh_id].cp)] = p_id_candidate2
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
    all_cp::Vector{Int},
    n_cp_hh::Int,
    model::ABM
    )

    # Check amount of cp of all hh, if insufficient replenish by randomly
    # drawing other suppliers based on weights
    for hh_id in all_hh

        if length(model[hh_id].cp) < n_cp_hh

            # Determine which bp are available
            n_add_cp = n_cp_hh - length(model[hh_id].cp)
            poss_cp = filter(p_id -> p_id ∉ model[hh_id].cp, all_cp)

            # Determine weights based on prices, sample and add
            weights = map(cp_id -> 1 / model[cp_id].p[end], poss_cp)
            add_cp = sample(poss_cp, Weights(weights), n_add_cp)
            model[hh_id].cp = vcat(model[hh_id].cp, add_cp)
        end
    end
end