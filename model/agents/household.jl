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
    skill::Float64                              # skill level of household

    # Income and wealth variables
    total_I::Float64 = L * skill                # total income from all factors
    labor_I::Float64 = L * skill                # income from labor
    capital_I::Float64 = 0.0                    # income from capital
    UB_I::Float64 = 0.0                         # income from unemployment benefits
    socben_I::Float64 = 0.0                     # income from social benefits (outside of UB)
    s::Float64 = 0.0                            # savings rate
    W::Float64 = 200                            # wealth or cash on hand

    # Consumption variables
    C::Float64 = 0.0                           # budget
    C_actual::Float64 = 0.0                    # actual spending on consumption goods
    cp::Vector{Int} = Vector{Int}()            # connected cp
    unsat_dem::Dict{Int, Float64} = Dict()     # unsatisfied demands
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
    for cp_id in hh.cp
        hh.unsat_dem[cp_id] = 0.0
    end
end


"""
Sets consumption budget based on current wealth level
"""
function set_consumption_budget_hh!(
    hh::Household,
    all_W::Vector{Float64},
    globalparam::GlobalParam,
    model::ABM
    )

    # Update average price level of bp and lp
    update_average_price_hh!(hh, model)

    # Compute consumption budget
    compute_consumption_budget_hh!(hh, maximum(all_W), minimum(all_W), globalparam.α_cp)

    # Reset actual spending to zero
    hh.C_actual = 0.0
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
    W_max::Float64,
    W_min::Float64,
    α_cp::Float64
    )

    # W_max = 800.
    # W_min = 50.
    if W_max - W_min < 1.
        W_max = W_min + 1.
    end

    # println(W_max," ", W_min)

    if hh.W > 0
        # Scale W between 0 and 100
        if hh.W == W_min
            hh.C = hh.W
        else
            W_scaled = 100 * (hh.W - W_min) / (W_max - W_min)
            # hh.C = (hh.W / W_scaled) * min(W_scaled^α_cp, W_scaled)

            hh.C = hh.P̄ * (hh.W / W_scaled) * min((W_scaled / hh.P̄)^α_cp, W_scaled / hh.P̄)
            hh.C = min(hh.C, hh.W)
        end
        # hh.C = min(hh.P̄[end] * (hh.W / hh.P̄[end])^α_cp, hh.W)

        # println("   ", hh.W, "  ", hh.C)

        # Compute savings rate
        hh.s = hh.total_I > 0 ? ((hh.total_I + hh.capital_I + hh.socben_I) - hh.C) / (hh.total_I + hh.capital_I + hh.socben_I) : 0.0 
    else
        hh.C = 0.0
        hh.s = 0.0
    end
end


"""
Household receives ordered cg and mutates balance
"""
function receive_ordered_goods_hh!(
    hh::Household,
    tot_sales::Float64,
    unsat_demand::Vector{Float64},
    hh_D::Vector{Float64},
    all_cp::Vector{Int},
    n_hh::Int
    )

    # Decrease wealth with total sold goods
    hh.W -= tot_sales
    hh.C_actual += tot_sales

    for cp_id in hh.cp
        i = cp_id - n_hh
        hh.unsat_dem[cp_id] = hh_D[i] > 0 ? unsat_demand[i] / hh_D[i] : 0.0
    end
end


"""
Updates satisfying wage wˢ and requested wage wʳ
"""
function update_sat_req_wage_hh!(
    hh::Household, 
    ϵ::Float64,
    ω::Float64, 
    w_min::Float64
    )

    # Update wˢ using adaptive rule
    # hh.wˢ = ω * hh.wˢ + (1 - ω) * (hh.employed ? hh.w[end] : w_min)
    if hh.employed
        hh.wˢ = max(w_min, hh.w[end])
    else
        # hh.wˢ = max(w_min, hh.wˢ * (1 - ϵ))
        hh.wˢ = max(w_min, hh.w[end] * (1 - ϵ * hh.T_unemp))
    end

    if hh.employed
        hh.wʳ = max(w_min, hh.w[end] * (1 + ϵ))
    else
        hh.wʳ = max(w_min, hh.wˢ)
    end
end


"""
Lets households get income, either from UB or wage
"""
function receiveincome_hh!(
    hh::Household, 
    amount::Float64;
    capgains::Bool=false,
    isUB::Bool=false,
    socben::Bool=false
    )

    # Add income to total income amount
    hh.total_I += amount
    hh.W += amount

    if capgains
        # Capital gains are added directly to the wealth, as they are added
        # at the end of the period
        hh.capital_I = amount

    elseif isUB
        hh.UB_I = amount
        hh.labor_I = 0.0
        hh.T_unemp += 1
    elseif socben
        # Transfer income can come from unemployment or other social benefits
        hh.socben_I = amount
    else
        # If employed, add to labor income, else add to transfer income
        hh.labor_I = amount
        hh.UB_I = 0.0
        shift_and_append!(hh.w, hh.w[end])
    end
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
    delete!(hh.unsat_dem, bankrupt_cp)
end


"""
Decides whether to switch to other cp
"""
function decide_switching_all_hh!(
    globalparam::GlobalParam,
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    all_p::Vector{Int},
    n_cp_hh::Int,
    model::ABM,
    to
    )

    for hh_id in all_hh
        # Check if demand was constrained and for chance of changing cp
        if length(model[hh_id].unsat_dem) > 0 && rand() < globalparam.ψ_Q

            # Pick a supplier to change, first set up weights inversely proportional
            # to supplied share of goods

            create_weights(hh::Household, cp_id::Int)::Float64 = hh.unsat_dem[cp_id] > 0 ? 1 / hh.unsat_dem[cp_id] : 0.0
            weights = map(cp_id -> create_weights(model[hh_id], cp_id), model[hh_id].cp)

            # Sample producer to replace
            p_id_replaced = sample(model[hh_id].cp, Weights(weights))[1]

            filter!(p_id -> p_id ≠ p_id_replaced, model[hh_id].cp)

            # Add new cp if list not already too long
            if (length(model[hh_id].cp) < n_cp_hh && length(model[hh_id].cp) < length(all_cp))
                
                # p_id_new = sample(setdiff(all_cp, model[hh_id].cp))
                p_id_new = sample(all_cp)

                count = 0
                while p_id_new ∈ model[hh_id].cp
                    p_id_new = sample(all_cp)

                    # Ugly way to avoid inf loop, will change later
                    count += 1
                    if count == 100
                        break
                    end
                end

                if p_id_new ∉ model[hh_id].cp
                    push!(model[hh_id].cp, p_id_new)
                    delete!(model[hh_id].unsat_dem, p_id_replaced)
                    model[hh_id].unsat_dem[p_id_new] = 0.0
                end
            end
        end

        # Check if household will look for a better price
        if rand() < globalparam.ψ_P


            if length(model[hh_id].cp) < length(all_cp)

                # Randomly select a supplier that may be replaced
                p_id_candidate1 = sample(model[hh_id].cp)

                # Randomly pick another candidate from same type and see if price is lower
                # Ugly sample to boost performance
                # p_id_candidate2 = sample(setdiff(all_cp, model[hh_id].cp))

                p_id_candidate2 = sample(all_cp)

                count = 0
                while p_id_candidate2 ∈ model[hh_id].cp
                    p_id_candidate2 = sample(all_cp)

                    # Ugly way to avoid inf loop, will change later
                    count += 1
                    if count == 100
                        break
                    end
                end
                
                # Replace supplier if price of other supplier is lower 
                if model[p_id_candidate2].p[end] < model[p_id_candidate1].p[end] && p_id_candidate2 ∉ model[hh_id].cp
                    # model[hh_id].cp[findfirst(x->x==p_id_candidate1, model[hh_id].cp)] = p_id_candidate2
                    filter!(p_id -> p_id ≠ p_id_candidate1, model[hh_id].cp)
                    push!(model[hh_id].cp, p_id_candidate2)

                    delete!(model[hh_id].unsat_dem, p_id_candidate1)
                    model[hh_id].unsat_dem[p_id_candidate2] = 0.0
                end
            else
                # If all producers known, throw out producer with highest price
                p_highest_price = maximum(cp_id -> model[cp_id].p[end], model[hh_id].cp)
                filter!(p_id -> p_id ≠ p_highest_price, model[hh_id].cp)
                delete!(model[hh_id].unsat_dem, p_highest_price)
            end
        end
    end
end


"""
Refills amount of bp and lp in amount is below minimum. Randomly draws suppliers
    inversely proportional to prices.
"""
function refillsuppliers_hh!(
    hh::Household,
    all_cp::Vector{Int},
    n_cp_hh::Int,
    model::ABM
    )

    if length(hh.cp) < n_cp_hh

        # Determine which bp are available
        n_add_cp = n_cp_hh - length(hh.cp)
        poss_cp = filter(p_id -> p_id ∉ hh.cp, all_cp)

        # Determine weights based on prices, sample and add
        weights = map(cp_id -> 1 / model[cp_id].p[end], poss_cp)
        add_cp = sample(poss_cp, Weights(weights), n_add_cp)
        hh.cp = vcat(hh.cp, add_cp)

        for cp_id in add_cp
            hh.unsat_dem[cp_id] = 0.0
        end
    end
end


"""
Samples wage levels of households from an empirical distribution.
"""
function sample_skills_hh(
    initparam::InitParam
    )::Vector{Float64}

    skills = []
    while length(skills) < initparam.n_hh
        s = rand(LogNormal(0.0, initparam.σ_hh_I)) * initparam.scale_hh_I
        if s < 2.5e5
            push!(skills, s)
        end
    end

    # Normalize skills
    skills = initparam.n_hh .* skills ./ sum(skills)

    return skills
end


"""
    reset_incomes_hh!(hh::Household)

Resets types incomes of household back to 0.0. Capital income is reset only before
    the new capital gains are sent.
"""
function resetincomes_hh!(
    hh::Household
    )

    hh.total_I = 0.0
    # hh.labor_I = 0.0
    # hh.capital_I = 0.0
    # hh.transfer_I = 0.0
end