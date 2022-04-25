function consumermarket_process!(
    all_hh::Vector{Int},
    all_cp::Vector{Int}, 
    all_W_hh::Vector{Float64},
    gov_struct::Government,
    global_param::GlobalParam,
    t::Int,
    model::ABM,
    to
    )

    # Make a dictionary with all cp and the inventory they own
    cp_inventories = Dict(cp_id => model[cp_id].p[end] * model[cp_id].N_goods for cp_id in all_cp)

    cp_with_inventory = Vector{Int}()
    for cp_id in all_cp
        if model[cp_id].N_goods > 0
            push!(cp_with_inventory, cp_id)
        end
    end

    # Let households set their budget and order their goods
    @inbounds for hh_id in all_hh
        # Set consumption budget and shares of good types and place orders
        set_consumption_budget_hh!(model[hh_id], all_W_hh, global_param, model)

        @timeit to "cp orders" cp_orders, cp_inventories, cp_with_inventory = place_orders_hh!(model[hh_id].cp, model[hh_id].C, 
                                                            cp_inventories, cp_with_inventory, global_param, model, to)

        for (cp_id, qp) in cp_orders
            push!(model[cp_id].order_queue, (hh_id, qp / model[cp_id].p[end]))
        end
    end

    # cp's handle order queue, send orders, households track which cp could not
    # supply the demand.
    for cp_id in all_cp
        @timeit to "send orders" send_ordered_goods_cp!(model[cp_id], t, model, to)
        reset_queue_cp!(model[cp_id])
    end
end


function cpmarket_matching(
    cm_dat::CMData,
    all_N::Vector{Float64},
    all_C::Vector{Float64},
    )

    # Normalize weights
    cm_dat.weights_sum .= sum(cm_dat.weights; dims=2)
    cm_dat.weights ./= cm_dat.weights_sum
    sold_out = []

    for _ in 1:3

        # Spread consumption budget according to weights
        cm_dat.C_spread .= cm_dat.all_C 
        cm_dat.C_spread .*= cm_dat.weights

        # Compute demand per cp and find which cp is sold out
        cm_dat.demand_per_cp .= max.(floor.(sum(cm_dat.C_spread; dims=1)[1,:], digits=6), 0)
        sold_out = findall(cm_dat.all_N .<= cm_dat.demand_per_cp)

        # Sell goods
        cm_dat.frac_sellable .= replace(min.(cm_dat.all_N ./ cm_dat.demand_per_cp, 1), -Inf=>0.0)
        cm_dat.C_spread .*= cm_dat.frac_sellable'
        cm_dat.transactions .+= cm_dat.C_spread

        # Update how much was sold per hh and cp
        cm_dat.sold_per_hh .+= floor.(sum(cm_dat.C_spread; dims=2)[:,1], digits=6)
        cm_dat.sold_per_cp .+= floor.(sum(cm_dat.C_spread; dims=1)[1,:], digits=6)

        # Update consumption budget C and inventory N
        cm_dat.all_C .-= cm_dat.sold_per_hh
        cm_dat.all_N .-= cm_dat.sold_per_cp

        # Set weights of sold-out producers to zero
        cm_dat.weights[:,sold_out] .= 0.0

        # Renormalize weights
        cm_dat.weights_sum .= sum(cm_dat.weights; dims=2)
        cm_dat.weights ./= cm_dat.weights_sum
        replace!(cm_dat.weights, NaN=>0.0)
    end

    # Assert if no hh overspent or cp sold too much
    # @assert all(cm_dat.sold_per_cp .<= all_N)
    # @assert all(cm_dat.sold_per_hh .<= all_C)

    return cm_dat.transactions
end


"""
Updates market share f of cp
"""
function update_marketshares_cm!(
    all_cp::Vector{Int}, 
    model::ABM
    )

    total_D = sum(cp_id -> model[cp_id].D[end], all_cp)

    for cp_id in all_cp
        cp = model[cp_id]
        f = cp.D[end] / total_D
        push!(model[cp_id].f, f)
    end
end


function fill_in_weights_cm!(
    cm_dat::CMData, 
    all_cp::Vector{Int}, 
    all_hh::Vector{Int},
    model::ABM
    )

    # First fill in a
end