function consumermarket_process!(
    all_hh::Vector{Int},
    all_cp::Vector{Int}, 
    all_bp::Vector{Int}, 
    all_lp::Vector{Int}, 
    all_W_hh::Vector{Float64},
    gov_struct::Government,
    global_param::GlobalParam,
    model::ABM
    )

    # Households set budget and send demand to queues of producers
    for hh_id in all_hh

        # Set consumption budget and shares of good types
        set_consumption_budget_hh!(
            model[hh_id],
            all_W_hh,
            global_param.c_L_max,
            global_param.a_σ,
            global_param.b_σ, 
            global_param.α_cp, 
            model
        )
        
        # Place orders at bp and lp
        place_orders_hh!(model[hh_id], model)
    end

    # cp's handle order queue, send orders, households track which cp could not
    # supply the demand.
    for cp_id in all_cp
        send_orders_cp!(model[cp_id], model)
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

    total_D = sum(map(cp_id -> model[cp_id].D[end], all_cp))

    for cp_id in all_cp
        cp = model[cp_id]
        f = cp.D[end] / total_D
        push!(model[cp_id].f, f)
    end
end