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

    total_wage = 0
    for hh_id in p.employees
        hh = model[hh_id]
        wage = hh.w[end] * hh.L
        total_wage += wage
        get_income_hh!(hh, wage)
    end
    
    p.curracc.TCL = total_wage
end


"""
Updates wage level for producer.
"""
function update_w̄_p!(
    p::AbstractAgent,
    model::ABM
    )

    if length(p.employees) > 0
        w̄ = mean(map(hh_id -> model[hh_id].w[end], p.employees))
        push!(p.w̄, w̄)
    end
end


"""
Updates market shares of cp firms
"""
function update_marketshare_cp!(
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    model::ABM
    )

    # Compute total market size of bp and lp markets
    bp_market = sum(bp_id -> model[bp_id].D[end], all_bp)
    lp_market = sum(lp_id -> model[lp_id].D[end], all_lp)

    # Update market share f for all bp
    for bp_id in all_bp
        if bp_market == 0
            f = 1/length(all_bp)
        else
            f = model[bp_id].D[end] / bp_market
        end
        push!(model[bp_id].f, f)
    end

    # Update market share f for all lp
    for lp_id in all_lp
        if lp_market == 0
            f = 1/length(all_lp)
        else
            f = model[lp_id].D[end] / lp_market
        end
        push!(model[lp_id].f, f)
    end
end


"""
Borrows amount of debt, permutates balance.
"""
function borrow_funds_p!(
    p::AbstractAgent,
    amount::Float64
    )

    # TODO permutate bank
    # p.balance.Deb += amount
    # p.balance.NW += amount

    # TODO: make this depend on amount of payback periods
    p.Deb_installments[1] += amount/3
    p.Deb_installments[2] += amount/3
    p.Deb_installments[3] += amount/3

end


"""
Pays back amount of debt, permutates balance.
"""
function payback_debt_p!(
    p::AbstractAgent,
    )

    # TODO permutate bank
    # p.balance.Deb -= amount
    # p.balance.NW -= amount

    to_be_paid = p.Deb_installments[1]
    p.curracc.rep_Deb += to_be_paid

    # TODO: make this depend on amount of payback periods
    p.Deb_installments[1] = p.Deb_installments[2]
    p.Deb_installments[2] = p.Deb_installments[3]
    p.Deb_installments[3] = 0.0
end

