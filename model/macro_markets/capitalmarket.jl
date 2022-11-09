mutable struct kmdata

    choices::Matrix{Float64}
    orders::Matrix{Float64}

    kp_capacity::Vector{Float64}

end



function capitalmarket_process!(
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    model::ABM
)

    # First gather all cp ids of cp that want to buy capital goods
    cp_ids = gather_expansionary_cp(all_cp, model)

    # Fill in first choices
    # TODO
    # TODO: have choices struct

    # Gather all kp ids that have production capacity (will be practically all at start)
    kp_ids = collect(1:n_kp)

    # Gather all production capacities of kp
    kp_cap = [5, 4, 4]

    # Empty orders data struct 
    # TODO: have orders struct
    orders = zeros(Int64, (3, 4))

    for roundnr in 1:1

        for choices in [firstchoices, secondchoices, thirdchoices]

            # Gather total faced demand
            faced_demand = map(kp_id -> sum(choices[kp_id, :]), kp_ids)

            println(faced_demand)

            # Sort kp ids based on demand faced in market
            # kp_ids = kp_ids[sortperm(faced_demand, rev=true)]

            # println(kp_ids)

            # Sort choices on size
            rankedchoices = sortslices(deepcopy(choices), by=x->sum(x), dims=1, rev=true)

            # Loop ONLY over kp ids that still have products
            for kp_id in kp_ids[sortperm(faced_demand, rev=true)]

                println(kp_id)
                # capacity = kp_cap[kp_i]

                # Loop ONLY over cp id that want to buy capital goods
                for cp_id in shuffle(1:n_cp)

                    if kp_cap[kp_id] == 0

                        # Take out kp id of kp id list
                        filter!(x -> x ≠ kp_id, kp_ids)

                        break
                    end

                    desired = max(choices[kp_id, cp_id] - sum(orders[:, cp_id]), 0)
                    order = min(kp_cap[kp_id], desired)

                    kp_cap[kp_id] -= order
                    orders[kp_id, cp_id] += order

                    # If full demand satisfied, delete cp id
                    # TODO
                end
            end
            display(orders)
        end

        # Recompute next choices based on now placed orders
        # TODO
    end

end


function gather_expansionary_cp(
    all_cp::Vector{Int64},
    model::ABM
)
    cp_expansionary = Int64[]

    for cp_id in all_cp
        if model[cp_id].Iᵈ > 0
            push!(cp_expansionary, cp_id)
        end
    end

    return cp_expansionary
end