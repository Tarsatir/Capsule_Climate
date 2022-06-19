function consumermarket_process!(
    all_hh::Vector{Int},
    all_cp::Vector{Int}, 
    government::Government,
    globalparam::GlobalParam,
    cm_dat::CMData,
    t::Int,
    model::ABM,
    to
    )

    # Households set consumption budget
    @timeit to "set budget" @inbounds for hh_id in all_hh
        set_consumption_budget_hh!(model[hh_id], globalparam, model)
    end

    # Reset cm data 
    @timeit to "reset matrices" reset_matrices_cp!(cm_dat, all_hh, all_cp, model)

    # Market clearing process
    @timeit to "market clearing" cpmarket_matching_cp!(cm_dat)

    @timeit to "process transac" process_transactions_cm!(all_hh, all_cp, cm_dat, model, to)

end


function process_transactions_cm!(
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    cm_dat::CMData, 
    model::ABM,
    to
    )

    unsat_demand = zeros(Float64, length(all_cp))
    hh_D = zeros(Float64, length(all_cp))

    # Process transactions for hh
    for (i, hh_id) in enumerate(minimum(all_hh):maximum(all_hh))

        hh_D .= @view cm_dat.true_D[i,:]

        # Compute unsatisfied demand
        unsat_demand .= @view cm_dat.true_D[i,:]
        unsat_demand .-= @view cm_dat.transactions[i,:]
        unsat_demand .= max.(unsat_demand, 0.0)

        receive_ordered_goods_hh!(
            model[hh_id], 
            sum(@view cm_dat.transactions[i,:]), 
            unsat_demand, 
            hh_D, 
            all_cp,
            length(all_hh)
        )
    end

    unsat_demand = zeros(Float64, length(all_hh))

    # Process transations for cp
    for (i, cp_id) in enumerate(minimum(all_cp):maximum(all_cp))

        # Compute unsat_demand
        unsat_demand = @view cm_dat.true_D[:,i]
        unsat_demand .-= @view cm_dat.transactions[:,i]
        unsat_demand .= max.(unsat_demand, 0.0)

        # TODO decide whether it makes sense if producers know unsatisfied demand

        model[cp_id].curracc.S = sum(@view cm_dat.transactions[:,i])
        N_goods_sold = model[cp_id].curracc.S / model[cp_id].p[end]
        shift_and_append!(model[cp_id].D, N_goods_sold)
        shift_and_append!(model[cp_id].Dᵁ, sum(unsat_demand))
        model[cp_id].N_goods = abs(model[cp_id].N_goods - N_goods_sold) > 1e-1 ? model[cp_id].N_goods - N_goods_sold : 0.0
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
    sold_out = Int64[]

    # println("   total N: $(sum(cm_dat.all_N))")
    # println("   total C: $(sum(cm_dat.all_C))")
    # println(cm_dat.all_N)

    # display(cm_dat.all_C)

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
        cm_dat.all_N[cm_dat.all_N .<= 1e-1] .= 0.0

        # Set weights of sold-out producers to zero
        cm_dat.weights[:,sold_out] .= 0.0

        # Renormalize weights
        sum!(cm_dat.weights_sum, cm_dat.weights)

        cm_dat.weights ./= cm_dat.weights_sum
        replace!(cm_dat.weights, NaN=>0.0)

        cm_dat.sold_per_hh_round .= 0.0
        cm_dat.sold_per_cp_round .= 0.0
    end

    # display(cm_dat.all_C)

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