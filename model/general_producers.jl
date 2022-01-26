function fire_excess_workers_p!(p)
    n_to_be_fired = abs(floor(Int, p.ΔLᵈ / 100))

    # TODO: find a more sophisticated way to select who is fired
    fired_workers = sample(p.Emp, n_to_be_fired, replace=false)

    # remove employees from labor stock
    p.L -= sum(map(hh -> hh.L, fired_workers))
    filter!(e -> e ∉ fired_workers, p.Emp)

    return fired_workers

end


function hire_worker_p!(p, l)

    # update labor stock and desired labor
    push!(p.Emp, l)
    p.L += l.L
    p.ΔLᵈ -= l.L

end


function pay_workers_p!(p)

    # loop over workers and pay out wage
    for hh in p.Emp
        push!(hh.I, hh.w[end] * hh.L)
    end

end