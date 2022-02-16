function fire_excess_workers_p!(
    p::AbstractAgent, 
    model::ABM
    )::Vector{Int}

    n_to_be_fired = floor(Int, abs(p.ΔLᵈ / 100))

    if n_to_be_fired >= length(p.employees)
        fired_workers = p.employees
    else
        # TODO: find a more sophisticated way to select who is fired
        fired_workers = sample(p.employees, n_to_be_fired, replace=false)
    end

    # println(fired_workers, length(fired_workers))

    # Remove employees from labor stock
    if length(fired_workers) > 0

        for hh_id in fired_workers
            remove_worker_p!(p, model[hh_id])
        end

        return fired_workers
    else
        return Vector{Int}()
    end

end


"""
Hires workers, adds worker id to worker array, updates labor stock.
"""
function hire_worker_p!(
    p::AbstractAgent, 
    hh::Household
    )

    # update labor stock and desired labor
    push!(p.employees, hh.id)
    p.L += hh.L
    p.ΔLᵈ -= hh.L
end


"""
Removes worker from firm
"""
function remove_worker_p!(
    p::AbstractAgent,
    hh::Household
    )

    p.L -= hh.L
    p.employees = filter(hh_id -> hh_id ≠ hh.id, p.employees)
end


"""
Loop over workers and pay out wage.
"""
function pay_workers_p!(
    p::AbstractAgent,
    model::ABM
    )
    for hh_id in p.employees
        hh = model[hh_id]
        total_wage = hh.w[end] * hh.L
        # if isnan(total_wage)
        #     println(hh.w[end], " ", hh.L, " ", total_wage)
        # end
        get_income_hh!(hh, total_wage)
    end
end


"""
Updates wage level for producer.
"""
function update_wage_level_p!(
    p::AbstractAgent,
    model::ABM
    )
    if length(p.employees) > 0
        w̄ = mean(map(hh_id -> model[hh_id].w[end], p.employees))
        push!(p.w̄, w̄)
    end
end