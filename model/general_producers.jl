function fire_excess_workers_p!(p :: AbstractAgent)

    n_to_be_fired = floor(Int, abs(p.ΔLᵈ / 100))

    # TODO: find a more sophisticated way to select who is fired
    if n_to_be_fired >= length(p.Emp)
        fired_workers = p.Emp
    else
        fired_workers = sample(p.Emp, n_to_be_fired, replace=false)
    end

    # remove employees from labor stock
    if length(fired_workers) > 0
        p.L -= sum(map(hh -> hh.L, fired_workers))
        filter!(e -> e ∉ fired_workers, p.Emp)

        p.P_FE = length(fired_workers) / length(p.Emp)

        return fired_workers
    else
        return []
    end

end


function hire_worker_p!(p :: AbstractAgent, l :: AbstractAgent)

    # update labor stock and desired labor
    push!(p.Emp, l)
    p.L += l.L
    p.ΔLᵈ -= l.L
    # println(p.ΔLᵈ)

end


function pay_workers_p!(p :: AbstractAgent)

    # loop over workers and pay out wage
    for hh in p.Emp
        
        total_wage = hh.w[end] * hh.L
        get_income_hh!(hh, total_wage)

    end

end


function update_wage_level_p!(p :: AbstractAgent)
    if length(p.Emp) > 0
        w̄ = mean(map(l -> l.w[end], p.Emp))
        push!(p.w, w̄)
    end
end