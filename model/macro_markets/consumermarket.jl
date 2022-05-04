function consumermarket_process!(
    all_hh::Vector{Int},
    all_cp::Vector{Int}, 
    all_W_hh::Vector{Float64},
    gov_struct::Government,
    global_param::GlobalParam,
    cm_dat::CMData,
    t::Int,
    model::ABM,
    to
    )

    # Households set consumption budget
    @inbounds for hh_id in all_hh
        @timeit to "set budget" set_consumption_budget_hh!(model[hh_id], all_W_hh, global_param, model)
    end

    # Reset cm data 
    @timeit to "reset matrices" reset_matrices_cp!(cm_dat, all_hh, all_cp, model)

    # for (cp_i, N) in enumerate(cm_dat.all_N)
    #     if N == 0
    #         cp = model[all_cp[cp_i]]
    #         println("age: $(cp.age), L: $(cp.L), ΔL: $(cp.ΔLᵈ), π_LP: $(cp.π_LP), mach: $(cp.n_machines), Q: $(cp.Q[end]), Qˢ: $(cp.Qˢ), Nᵈ: $(cp.Nᵈ), Dᵉ: $(cp.Dᵉ), D: $(cp.D[end])")
    #     end
    # end

    # Market clearing process
    @timeit to "market clearing" cpmarket_matching_cp!(cm_dat)

    @timeit to "process transac" process_transactions_cm!(all_hh, all_cp, cm_dat, model, to)

    # TODO: so smt with actual demand

    # @timeit to "write demand" write_demand_hh_to_cp!(all_cp, all_hh, transactions, model)

    # @timeit to "send orders" send_ordered_goods_all_cp!(all_cp, t, model, to)

    # # Make a dictionary with all cp and the inventory they own
    # cp_inventories = Dict(cp_id => model[cp_id].p[end] * model[cp_id].N_goods for cp_id in all_cp)

    # cp_with_inventory = Vector{Int}()
    # for cp_id in all_cp
    #     if model[cp_id].N_goods > 0
    #         push!(cp_with_inventory, cp_id)
    #     end
    # end

    # # Let households set their budget and order their goods
    # @inbounds for hh_id in all_hh
    #     # Set consumption budget and shares of good types and place orders
    #     set_consumption_budget_hh!(model[hh_id], all_W_hh, global_param, model)

        # @timeit to "cp orders" cp_orders, cp_inventories, cp_with_inventory = place_orders_hh!(model[hh_id].cp, model[hh_id].C, 
    #                                                         cp_inventories, cp_with_inventory, global_param, model, to)

    #     for (cp_id, qp) in cp_orders
    #         push!(model[cp_id].order_queue, (hh_id, qp / model[cp_id].p[end]))
    #     end
    # end

    # cp's handle order queue, send orders, households track which cp could not
    # supply the demand.
    # for cp_id in all_cp
        # @timeit to "send orders" send_ordered_goods_cp!(model[cp_id], t, model, to)
    #     reset_queue_cp!(model[cp_id])
    # end
end


function process_transactions_cm!(
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    cm_dat::CMData, 
    model::ABM,
    to
    )

    unsat_demand = zeros(length(all_cp))
    hh_D = zeros(length(all_cp))

    # Process transactions for hh
    for (i, hh_id) in enumerate(all_hh)

        # hh_transac .= cm_dat.transactions[i,:]
        hh_D .= cm_dat.true_D[i,:]

        # Compute unsatisfied demand
        @timeit to "c1" unsat_demand = hh_D
        @timeit to "c2" unsat_demand .-= @view cm_dat.transactions[i,:]
        unsat_demand .= max.(unsat_demand, 0.0)

        @timeit to "receive goods" receive_ordered_goods_hh!(
                                        model[hh_id], 
                                        sum(@view cm_dat.transactions[i,:]), 
                                        unsat_demand, 
                                        hh_D, 
                                        all_cp
                                    )
    end

    unsat_demand = zeros(length(all_hh))
    # cp_transac = zeros(length(all_hh))
    # cp_D = zeros(length(all_hh))

    # Process transations for cp
    for (i, cp_id) in enumerate(all_cp)

        # Compute unsat_demand
        @timeit to "c3" unsat_demand = @view cm_dat.true_D[:,i]
        @timeit to "c4" unsat_demand .-= @view cm_dat.transactions[:,i]
        unsat_demand .= max.(unsat_demand, 0.0)

        # TODO decide whether it makes sense if producers know unsatisfied demand

        @timeit to "c5" model[cp_id].curracc.S = sum(@view cm_dat.transactions[:,i])
        N_goods_sold = model[cp_id].curracc.S / model[cp_id].p[end]
        shift_and_append!(model[cp_id].D, N_goods_sold)
        shift_and_append!(model[cp_id].Dᵁ, sum(unsat_demand))
        model[cp_id].N_goods = abs(model[cp_id].N_goods - N_goods_sold) < 1e-1 ? model[cp_id].N_goods - N_goods_sold : 0.0
    end

end


"""
Market matching process for consumer market.
"""
function cpmarket_matching_cp!(
    cm_dat::CMData
    )

    # Normalize weights
    sum!(cm_dat.weights_sum, cm_dat.weights)
    cm_dat.weights ./= cm_dat.weights_sum
    sold_out = []

    # println("total N: $(sum(cm_dat.all_N))")
    # println(cm_dat.all_N)

    for i in 1:3

        # Spread consumption budget according to weights (no allocs)
        cm_dat.C_spread .= cm_dat.all_C 
        cm_dat.C_spread .*= cm_dat.weights

        # Compute the first choice of hh for the products (2 allocs)
        if i == 1
            cm_dat.true_D .= cm_dat.C_spread
        end

        # Compute demand per cp and find which cp is sold out
        sum!(cm_dat.demand_per_cp, cm_dat.C_spread')
        cm_dat.demand_per_cp .= max.(floor.(cm_dat.demand_per_cp, digits=8), 0.0)

        sold_out = findall(cm_dat.all_N .<= cm_dat.demand_per_cp)

        # Sell goods
        cm_dat.frac_sellable .= min.(cm_dat.all_N ./ cm_dat.demand_per_cp, 1.0)
        replace!(cm_dat.frac_sellable, -Inf=>0.0, NaN=>0.0)

        cm_dat.C_spread .*= cm_dat.frac_sellable'
        cm_dat.transactions .+= cm_dat.C_spread

        # Update how much was sold per hh and cp
        sum!(cm_dat.sold_per_hh_round, cm_dat.C_spread)
        # cm_dat.sold_per_hh .+= floor.(cm_dat.sold_per_hh_round, digits=8)
        
        sum!(cm_dat.sold_per_cp_round, cm_dat.C_spread')
        # cm_dat.sold_per_cp .+= floor.(cm_dat.sold_per_cp_round, digits=6)

        # Update consumption budget C and inventory N
        cm_dat.all_C .-= cm_dat.sold_per_hh_round
        cm_dat.all_C[cm_dat.all_C .<= 1e-1] .= 0.0

        cm_dat.all_N .-= cm_dat.sold_per_cp_round
        cm_dat.all_C[cm_dat.all_C .<= 1e-1] .= 0.0

        # Set weights of sold-out producers to zero
        cm_dat.weights[:,sold_out] .= 0.0

        # Renormalize weights
        sum!(cm_dat.weights_sum, cm_dat.weights)

        cm_dat.weights ./= cm_dat.weights_sum
        replace!(cm_dat.weights, NaN=>0.0)

        cm_dat.sold_per_hh_round .= 0.0
        cm_dat.sold_per_cp_round .= 0.0
    end

    # println("Total sales: $(sum(cm_dat.transactions))")

    # return cm_dat.transactions
end


"""
Pushes demanded goods of hh to cp.
"""
function write_demand_hh_to_cp!(
    all_cp::Vector{Int}, 
    all_hh::Vector{Int}, 
    transactions::Matrix{Float64}, 
    model::ABM
    )

    for (cp_id, order_col) in zip(all_cp, eachcol(transactions))
        for (hh_id, qp) in zip(all_hh, order_col)
            push!(model[cp_id].order_queue, (hh_id, qp / model[cp_id].p[end]))
        end
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


function fill_in_weights_cm!(
    cm_dat::CMData, 
    all_cp::Vector{Int}, 
    all_hh::Vector{Int},
    model::ABM
    )

    # First fill in a
end