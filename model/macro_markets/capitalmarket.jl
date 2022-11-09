mutable struct KMData

    choices::Matrix{Int64}                  # Matrix containing choices of cp
    orders::Matrix{Int64}                   # Matrix containing orders of cp

    kp_capacity::Vector{Int64}              # Vector of production capacity kp

end


function capitalmarket_process!(
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    model::ABM
)

    # First gather all cp ids of cp that want to buy capital goods
    cp_ids = gather_expansionary_cp(all_cp, model)

    # Gather all kp ids and production capacity
    gather_capacities_kp!(all_kp, model)

    # println(cp_ids)
    # println(kp_capacity)

    # Fill in first choices
    fill_in_choices!(all_cp, all_kp, model)

    # Empty orders data struct 
    # TODO

    # for roundnr in 1:1

    #     for choices in [firstchoices, secondchoices, thirdchoices]

    #         # Gather total faced demand
    #         faced_demand = map(kp_id -> sum(choices[kp_id, :]), kp_ids)

    #         println(faced_demand)

    #         # Sort kp ids based on demand faced in market
    #         # kp_ids = kp_ids[sortperm(faced_demand, rev=true)]

    #         # println(kp_ids)

    #         # Sort choices on size
    #         rankedchoices = sortslices(deepcopy(choices), by=x->sum(x), dims=1, rev=true)

    #         # Loop ONLY over kp ids that still have products
    #         for kp_id in kp_ids[sortperm(faced_demand, rev=true)]

    #             println(kp_id)
    #             # capacity = kp_cap[kp_i]

    #             # Loop ONLY over cp id that want to buy capital goods
    #             for cp_id in shuffle(1:n_cp)

    #                 if kp_cap[kp_id] == 0

    #                     # Take out kp id of kp id list
    #                     filter!(x -> x ≠ kp_id, kp_ids)

    #                     break
    #                 end

    #                 desired = max(choices[kp_id, cp_id] - sum(orders[:, cp_id]), 0)
    #                 order = min(kp_cap[kp_id], desired)

    #                 kp_cap[kp_id] -= order
    #                 orders[kp_id, cp_id] += order

    #                 # If full demand satisfied, delete cp id
    #                 # TODO
    #             end
    #         end
    #         display(orders)
    #     end

    #     # Recompute next choices based on now placed orders
    #     # TODO
    # end

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


function fill_in_choices!(
    all_cp::Vector{Int64}, 
    all_kp::Vector{Int64}, 
    model::ABM
)

    # TODO

end


"""
Gathers production capacities of all kp
"""
function gather_capacities_kp!(
    all_kp::Vector{Int64}, 
    model::ABM
)

    model.kmdata.kp_capacity .= map(kp_id -> model[kp_id].prod_cap, all_kp)
end