function consumermarket_process!(
    all_hh::Vector{Int},
    all_cp::Vector{Int}, 
    all_bp::Vector{Int}, 
    all_lp::Vector{Int}, 
    all_W_hh::Vector{Float64},
    gov_struct::Government,
    global_param::GlobalParam,
    t::Int,
    model::ABM,
    to
    )

    # Make a dictionary with all cp and the inventory they own
    cp_inventories = Dict(cp_id => model[cp_id].p[end] * model[cp_id].N_goods for cp_id in all_cp)

    bp_with_inventory = Vector{Int}()
    for bp_id in all_bp
        if model[bp_id].N_goods > 0
            push!(bp_with_inventory, bp_id)
        end
    end

    lp_with_inventory = Vector{Int}()
    for lp_id in all_lp
        if model[lp_id].N_goods > 0
            push!(lp_with_inventory, lp_id)
        end
    end


    # Let households set their budget and order their goods
    @inbounds for hh_id in all_hh
        # Set consumption budget and shares of good types and place orders
        set_consumption_budget_hh!(model[hh_id], all_W_hh, global_param, model)

        # Divide budget by chosen category
        C_b = model[hh_id].C[end] * (1 - model[hh_id].c_L)
        C_l = model[hh_id].C[end] - C_b

        @timeit to "bp orders" bp_orders, cp_inventories, bp_with_inventory = place_orders_hh!(model[hh_id].bp, C_b, cp_inventories, 
                                                         bp_with_inventory, global_param, model, to)
        @timeit to "lp orders" lp_orders, cp_inventories, lp_with_inventory = place_orders_hh!(model[hh_id].lp, C_l, cp_inventories, 
                                                         lp_with_inventory, global_param, model, to)
        # Send order to queues of bp and lp
        for (cp_id, qp) in merge(bp_orders, lp_orders)
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