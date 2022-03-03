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
    # p.ΔLᵈ -= hh.L
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
Adds borrowed amount as an incoming cashflow to current account.
"""
function borrow_funds_p!(
    p::AbstractAgent,
    amount::Float64,
    b=3::Int64
    )

    # Add debt as repayments for coming periods
    for t in 2:b+1
        p.debt_installments[t] += amount / b
    end

    # Add received funds as incoming cashflow
    p.curracc.add_debt += amount
end


"""
Adds to-be-repaid amount as an outgoing cashflow to current account.
"""
function payback_debt_p!(
    p::AbstractAgent,
    b=3::Int64
    )

    # Add repaid debt as outgoing cashflow
    p.curracc.rep_debt = p.debt_installments[1]

    # Shift remaining debt amounts.
    for t in 1:b
        p.debt_installments[t] = p.debt_installments[t + 1]
    end
    p.debt_installments[b+1] = 0.0

    p.balance.debt -= p.curracc.rep_debt
end

