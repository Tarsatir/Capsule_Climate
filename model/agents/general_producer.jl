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
    # t::Int,
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
    # t::Int,
    model::ABM
    )

    if length(p.employees) > 0
        # p.w̄[end] = mean(hh_id -> model[hh_id].w[end], p.employees)
        shift_and_append!(p.w̄, mean(hh_id -> model[hh_id].w[end], p.employees))
    else
        # p.w̄[end] = p.w̄[t-1]
        shift_and_append!(p.w̄, p.w̄[end])
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
            f = 1 / length(all_bp)
        else
            f = model[bp_id].D[end] / bp_market
        end
        shift_and_append!(model[bp_id].f, f)
    end

    # Update market share f for all lp
    for lp_id in all_lp
        if lp_market == 0
            f = 1/length(all_lp)
        else
            f = model[lp_id].D[end] / lp_market
        end
        shift_and_append!(model[lp_id].f, f)
    end
end


"""
Adds borrowed amount as an incoming cashflow to current account.
"""
function borrow_funds_p!(
    p::AbstractAgent,
    amount::Float64,
    b::Int64
    )

    # Add debt as repayments for coming periods
    for i in 2:b+1
        p.debt_installments[i] += amount / b
    end

    # Add received funds as incoming cashflow
    p.curracc.add_debt += amount
    p.balance.debt = sum(p.debt_installments)
end


"""
Adds to-be-repaid amount as an outgoing cashflow to current account.
"""
function payback_debt_p!(
    p::AbstractAgent,
    b::Int64
    )

    # println("debt: $(p.balance.debt), debt installments: $(p.debt_installments)")

    # Add repaid debt as outgoing cashflow
    p.curracc.rep_debt += p.debt_installments[1]

    # Shift remaining debt amounts.
    for i in 1:b
        p.debt_installments[i] = p.debt_installments[i+1]
    end
    p.debt_installments[b+1] = 0.0

    p.balance.debt = sum(p.debt_installments)
end


"""
Checks whether producers are bankrupt and should be replaced, returns vectors
    containing ids of to-be-replaced producers by type
"""
function check_bankrupty_all_p!(
    all_p::Vector{Int},
    all_kp::Vector{Int},
    global_param::GlobalParam, 
    model::ABM
    )::Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}}

    bankrupt_bp = Vector{Int}()
    bankrupt_lp = Vector{Int}()
    bankrupt_kp = Vector{Int}()
    bankrupt_kp_i = Vector{Int}()

    # return bankrupt_bp, bankrupt_lp, bankrupt_kp, bankrupt_kp_i

    bp_counter = 0
    lp_counter = 0
    kp_counter = 0

    for p_id in all_p
        if check_if_bankrupt_p!(model[p_id], global_param.t_wait)
            if typeof(model[p_id]) == ConsumerGoodProducer
                if model[p_id].type_good == "Basic"
                    bp_counter += 1
                    push!(bankrupt_bp, p_id)
                else
                    lp_counter += 1
                    push!(bankrupt_lp, p_id)
                end
            else
                kp_counter += 1
                if kp_counter != length(all_kp)
                    push!(bankrupt_kp, p_id)
                    push!(bankrupt_kp_i, model[p_id].kp_i)
                end
            end
        end
    end

    # println("Bankrupties, kp: $kp_counter, bp: $bp_counter, lp: $lp_counter")

    return bankrupt_bp, bankrupt_lp, bankrupt_kp, bankrupt_kp_i
end


"""
Kills all bankrupt producers.
    - for cp, removes bankrupt companies from households.
    - for all, fires all remaining workers.
    - for all, removes firm agent from model.
"""
function kill_all_bankrupt_p!(
    bankrupt_bp::Vector{Int},
    bankrupt_lp::Vector{Int},
    bankrupt_kp::Vector{Int},
    all_hh::Vector{Int},
    all_kp::Vector{Int},
    labormarket_struct,
    indexfund_struct,
    model::ABM
    )

    # Remove bankrupt cp ids from households
    for hh_id in all_hh
        remove_bankrupt_producers_hh!(model[hh_id], bankrupt_bp, bankrupt_lp)
    end

    # Remove bankrupt cp ids from kp historical clients
    for kp_id in all_kp
        remove_bankrupt_HC_kp!(model[kp_id], bankrupt_bp, bankrupt_lp)
    end

    # All employees are declared unemployed
    all_employees = Vector{Int}()
    for p_id in Iterators.flatten((bankrupt_bp, bankrupt_lp, bankrupt_kp))
        append!(all_employees, model[p_id].employees)
    end
    
    update_firedworker_lm!(labormarket_struct, all_employees)

    ages = []
    Qs = []

    # Fire all workers still remaining in the firm, remove firm
    total_unpaid_net_debt = 0.0
    for p_id in Iterators.flatten((bankrupt_bp, bankrupt_lp, bankrupt_kp))

        # Fire remaining workers
        for hh_id in model[p_id].employees
            set_unemployed_hh!(model[hh_id])
        end

        # TODO: TEMP SOLUTION, DESCRIBE IT WORKS
        # indexfund_struct.Assets += (model[p_id].balance.NW - model[p_id].balance.debt)

        total_unpaid_net_debt += (model[p_id].balance.debt - model[p_id].balance.NW)
        # total_unpaid_net_debt += model[p_id].balance.debt

        if p_id ∈ bankrupt_bp || p_id ∈ bankrupt_lp
            push!(ages, model[p_id].age)
            push!(Qs, model[p_id].Q[end])
        end

        # Remove firm agents from model
        kill_agent!(p_id, model)
    end

    deduct_unpaid_net_debts_if!(indexfund_struct, total_unpaid_net_debt)
    # println("All outflow: $(total_unpaid_net_debt)")
    # println("Age bankrupt cp: $(length(ages) > 0 ? mean(ages) : "yeet")")
    # println("Prod bankrupt cp: $(length(Qs) > 0 ? mean(Qs) : "yeet")")
end


"""
Checks if producers must be declared bankrupt
"""
function check_if_bankrupt_p!(
    p::AbstractAgent,
    t_wait::Int
    )::Bool

    # if p.age > t_wait && (p.f[end] <= 0.0001 || p.balance.EQ < 0)
    if p.age > t_wait && p.f[end] <= 0.0001
        return true
    end
    return false
end 


"""
Determines new productivity level, resulting from innovation process for both kp and ep.
"""
function update_techparam_p(
    TP_last::Float64, 
    global_param::GlobalParam;
    is_EF::Bool = false
    )::Float64

    κ_i = rand(Beta(global_param.α1, global_param.β1))
    κ_i = global_param.κ_lower + κ_i * (global_param.κ_upper - global_param.κ_lower)

    return is_EF ? TP_last*(1 - κ_i) : TP_last*(1 + κ_i)
end