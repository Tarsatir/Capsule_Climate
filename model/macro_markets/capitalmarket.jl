mutable struct KMData

    choices::Matrix{Int64}                  # Matrix containing choices of cp
    orders::Matrix{Int64}                   # Matrix containing orders of cp

    # kp_capacity::Vector{Int64}              # Vector of production capacity kp

end


function capitalmarket_process!(
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    government::Government,
    globalparam::GlobalParam,
    ep::EnergyProducer,
    t::Int64,
    model::ABM;
    total_n_rounds::Int64=1
)

    for cp_id in all_cp
        plan_replacement_cp!(model[cp_id], government, globalparam, ep, 1, t, model)
        plan_expansion_cp!(model[cp_id], globalparam, 1, model)
    end

    # First gather all cp ids of cp that want to buy capital goods
    cp_ids = gather_expansionary_cp(all_cp, model)

    # Gather all kp ids and production capacity
    kp_ids = gather_producing_kp(all_kp, model)
    kp_capacity = gather_capacities_kp!(kp_ids, model)

    # println(cp_ids)
    # println(kp_capacity)

    faced_demand = nothing

    # Empty orders data struct 
    model.kmdata.orders .= 0

    for roundnr in 1:total_n_rounds

        # Fill in first choices
        fill_in_cp_choices!(cp_ids, roundnr, model)

        if roundnr == 1
            # Gather total faced demand for all kp
            faced_demand = Dict(kp_id => sum(model.kmdata.choices[model[kp_id].kp_i, :]) for kp_id ∈ kp_ids)
            println("FACED DEMAND")
            println(faced_demand)
        end

        println("CAPACITY")
        println(kp_capacity)

        for kp_id in kp_ids

            kp_i = model[kp_id].kp_i
            price = model.kp_brochures[Symbol(kp_id)][:price]

            findall(x -> x != 0, model.kmdata.choices[kp_i,:])
            
            for cp_i in findall(x -> x != 0, model.kmdata.choices[kp_i,:])

                # Recover cp id
                cp_id = get(get(model.i_to_id, "cp", nothing), cp_i, nothing)
                
                if model.kmdata.choices[kp_i, cp_i] < kp_capacity[kp_id]

                    # If capacity larger than demand, place full order and decrease capacity
                    order = model.kmdata.choices[kp_i, cp_i]
                    model.kmdata.orders[kp_i, cp_i] = order
                    kp_capacity[kp_id] -= order
                    process_machine_order_cp!(model[cp_id], order, price)

                else
                    
                    # If capacity smaller than demand, place partial order and decrease capacity
                    order = kp_capacity[kp_id]
                    model.kmdata.orders[kp_i, cp_i] = order
                    kp_capacity[kp_id] = 0
                    process_machine_order_cp!(model[cp_id], order, price)

                    filter!(j -> j ≠ kp_id, kp_ids)
                    delete!(kp_capacity, kp_id)
                    break

                end
            end
        end

        # Recompute next choices based on now placed orders
        if roundnr != total_n_rounds
            for cp_id in cp_ids
                plan_replacement_cp!(model[cp_id], government, globalparam, ep, 1, t, model)
                plan_expansion_cp!(model[cp_id], globalparam, 1, model)
            end
        end
    end

    # display(model.kmdata.orders)

    process_machine_orders!(model, globalparam.freq_per_machine, faced_demand)
end


function gather_expansionary_cp(
    all_cp::Vector{Int64},
    model::ABM
)

    cp_expansionary = Int64[]

    for cp_id in all_cp
        if model[cp_id].n_mach_desired_RS + model[cp_id].n_mach_desired_EI > 0
            push!(cp_expansionary, cp_id)
        end
    end

    return cp_expansionary
end


function gather_producing_kp(
    all_kp::Vector{Int64}, 
    model::ABM
)

    kp_producing = Int64[]

    for kp_id in all_kp
        if model[kp_id].L > 0
            push!(kp_producing, kp_id)
        end
    end

    return kp_producing
end


"""
Gathers production capacities of all kp
"""
function gather_capacities_kp!(
    kp_ids::Vector{Int64}, 
    model::ABM
)

    return Dict(kp_id => model[kp_id].prod_cap for kp_id ∈ kp_ids)
end


function fill_in_cp_choices!(
    cp_ids::Vector{Int64}, 
    roundnr::Int64,
    model::ABM
)

    # Clear previous choices
    model.kmdata.choices .= 0

    # Loop over cp and fill in choices
    for cp_id in cp_ids

        cp = model[cp_id]
        kp_i = model[cp.kp_ids[roundnr]].kp_i

        # Fill in total number of desired machines of nth choice producer
        model.kmdata.choices[kp_i, cp.cp_i] = cp.n_mach_desired_EI + cp.n_mach_desired_RS 
    end
end


function process_machine_orders!(
    model::ABM,
    freq_per_machine::Int64,
    faced_demand::Dict
)

    # Get all placed orders
    for idx in findall(x -> x ≠ 0, model.kmdata.orders)
        kp_i = idx[1]
        cp_i = idx[2]

        kp_id = get(get(model.i_to_id, "kp", nothing), kp_i, nothing)
        cp_id = get(get(model.i_to_id, "cp", nothing), cp_i, nothing)

        order_machines_cp!(
                            model[cp_id], 
                            model[kp_id], 
                            model.kmdata.orders[kp_i, cp_i],
                            freq_per_machine
                          )

        total_orders_kp = sum(model.kmdata.orders[kp_i, :])

        # Set kp total faced demand
        model[kp_id].O_unmet = max(faced_demand[kp_id] - total_orders_kp, 0)
    end
end