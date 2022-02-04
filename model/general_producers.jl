function fire_excess_workers_p!(
    p::AbstractAgent, 
    model::ABM
    ) :: Vector{Int}

    n_to_be_fired = floor(Int, abs(p.ΔLᵈ / 100))

    if n_to_be_fired >= length(p.Emp)
        fired_workers = p.Emp
    else
        # TODO: find a more sophisticated way to select who is fired
        fired_workers = sample(p.Emp, n_to_be_fired, replace=false)
    end

    # remove employees from labor stock
    if length(fired_workers) > 0
        p.L -= sum(map(hh_id -> model[hh_id].L, fired_workers))
        filter!(hh_id -> hh_id ∉ fired_workers, p.Emp)

        p.P_FE = length(fired_workers) / length(p.Emp)

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
    push!(p.Emp, hh.id)
    p.L += hh.L
    p.ΔLᵈ -= hh.L

end


"""
Loop over workers and pay out wage.
"""
function pay_workers_p!(
    p::AbstractAgent,
    model::ABM
    )
    for hh_id in p.Emp
        hh = model[hh_id]
        total_wage = hh.w[end] * hh.L
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
    if length(p.Emp) > 0
        w̄ = mean(map(hh_id -> model[hh_id].w[end], p.Emp))
        push!(p.w, w̄)
    end
end