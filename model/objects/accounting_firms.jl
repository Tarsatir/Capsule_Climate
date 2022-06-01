"""
BALANCE SHEET
"""
@with_kw mutable struct Balance
    # Assets
    N::Float64 = 0.0            # Inventories
    K::Float64 = 0.0            # Capital
    NW::Float64 = 0.0           # Liquid assets

    # Liabilities
    debt::Float64 = 0.0         # Debt
    EQ::Float64 = 0.0           # Equity
end


"""
Closes balance by computing firm's equity.
If firm is insolvent, liquidate firm.
"""
function close_balance_all_p!(
    all_p::Vector{Int},
    global_param::GlobalParam,
    τᴾ::Float64,
    indexfund_struct,
    t::Int,
    model::ABM
    )

    total_dividends = 0.0

    all_max_NW = 0
    all_NW = 0

    for p_id in all_p

        model[p_id].balance.debt = max(0.0, model[p_id].balance.debt)

        # Compute interest payment
        update_interest_payment_p!(model[p_id], global_param.r)

        # Repay debts of period
        payback_debt_p!(model[p_id], global_param.b)

        # Update valuation of capital stock
        writeoffs = update_K_p!(model[p_id], global_param.η)

        # Compute profits
        compute_Π_p!(model[p_id]; writeoffs)

        # Update liquid assets NW
        update_NW_p!(model[p_id], τᴾ)

        # Update valuation of inventory. kp do not have inventory.
        if typeof(model[p_id]) == ConsumerGoodProducer
            model[p_id].balance.N = model[p_id].p[end] * model[p_id].N_goods
        end

        # Determine how much the firm can have as NW at most
        max_NW = global_param.max_NW_ratio * (model[p_id].curracc.TCL + model[p_id].curracc.TCI + 
                      model[p_id].curracc.int_debt + model[p_id].debt_installments[2])

        if typeof(model[p_id]) == ConsumerGoodProducer
            all_max_NW += max_NW
        end

        # If not enough liquid assets available, borrow additional funds.
        if model[p_id].balance.NW < 0
            # max_add_debt = max(model[p_id].curracc.S * Λ - model[p_id].balance.debt, 0)
            # add_debt = min(-model[p_id].balance.NW, max_add_debt)
            # borrow_funds_p!(model[p_id], add_debt)
            # model[p_id].balance.NW += add_debt
            borrow_funds_p!(model[p_id], -model[p_id].balance.NW, global_param.b)
            model[p_id].balance.NW = 0
        elseif (!check_if_bankrupt_p!(model[p_id],  global_param.t_wait) 
                && (model[p_id].balance.NW > max_NW) && (t > global_param.t_wait))
        # elseif (model[p_id].balance.NW > max_NW && t > 1)
            # indexfund_struct.Assets += (model[p_id].balance.NW - max_NW)
            total_dividends += model[p_id].balance.NW - max_NW
            model[p_id].balance.NW = max_NW
        end

        # Compute Equity
        tot_assets = model[p_id].balance.N + model[p_id].balance.K + model[p_id].balance.NW
        model[p_id].balance.EQ = tot_assets - model[p_id].balance.debt

        # If NW is negative, maximum debt is reached, and EQ is set to
        # a negative value so the firm is declared bankrupt
        # if model[p_id].balance.debt > global_param.Λ * model[p_id].curracc.S
        #     model[p_id].balance.EQ = -1.0
        #     model[p_id].f[end] = 0.0
        # end

        if typeof(model[p_id]) == ConsumerGoodProducer
            all_NW += model[p_id].balance.NW
        end
    end

    # println("avg max NW: ", all_max_NW / 200, ", all max NW: ", all_max_NW)
    # println("avg NW stock: ", all_NW / 200, ", all NW: ", all_NW)

    # println("total dividends: ", total_dividends)
    receive_dividends_if!(indexfund_struct, total_dividends)
end


"""
Computes current value of capital stock and writeoffs in period
"""
function update_K_p!(
    p::AbstractAgent,
    η::Int
    )::Float64

    # Only cp have capital stock
    if typeof(p) == ConsumerGoodProducer
        # Loop over machines, add current value of machines
        K = 0.0
        writeoffs = 0.0
        for machine in p.Ξ
            newval = machine.p * machine.freq
            K += newval
            # K += max((η - machine.age) / η, 0) * newval
            # if machine.age <= η
            #     writeoffs += (1 / η) * newval
            # end
        end

        p.balance.K = K

        return writeoffs
    else
        return 0.0
    end
end


"""
Updates NW of firm based on last period's spending.
"""
function update_NW_p!(
    p::AbstractAgent,
    τᴾ::Float64
    )

    NW = p.balance.NW
    S = p.curracc.S
    TCL = p.curracc.TCL
    TCI = p.curracc.TCI
    TCE = p.curracc.TCE
    rep_debt = p.curracc.rep_debt
    add_debt = p.curracc.add_debt
    int_debt = p.curracc.int_debt
    rev_dep = p.curracc.rev_dep
    profit_tax = max(τᴾ * p.Π[end], 0)

    # if typeof(p) == ConsumerGoodProducer
    #     println("NW: $NW, S: $S, TCL: $TCL, TCI: $TCI, rd: $rep_debt, ad: $add_debt, id: $int_debt, pt: $profit_tax")
    # end

    new_NW = NW + S + add_debt + rev_dep - TCL - TCI - TCE - rep_debt - int_debt - profit_tax

    if typeof(p) == ConsumerGoodProducer
        p.NW_growth = (new_NW - p.balance.NW) / p.balance.NW
    end
    
    p.balance.NW = new_NW
end


"""
Computes the amount of interest the firm has to pay over its debts
"""
function update_interest_payment_p!(
    p::AbstractAgent,
    r::Float64
    )

    p.curracc.int_debt = r * p.balance.debt
end


"""
Computes profit of cp in previous time period
"""
function compute_Π_p!(
    p::AbstractAgent;
    writeoffs=0.0::Float64
    )
 
    Π = p.curracc.S + p.curracc.rev_dep - p.curracc.TCL - p.curracc.TCE - p.curracc.int_debt - writeoffs
    # if typeof(p) == CapitalGoodProducer
    #     println("profit ", Π, " S ", p.curracc.S, "  TCL ", p.curracc.TCL, " TCI ", p.curracc.TCI, " int ", p.curracc.int_debt, " p ", p.p[end])
    #     println("labor ", p.ΔLᵈ, " O ", p.O, " O/B ", p.O/p.B_LP[end], " RD ", p.RD / p.w̄[end], " L ", p.L, " w ", p.w̄[end])
    # end
    # Π = p.curracc.S + p.curracc.rev_dep - p.curracc.TCL - p.curracc.int_debt - writeoffs

    # push!(p.Π, Π)
    shift_and_append!(p.Π, Π)
end


"""
CURRENT ACCOUNT
"""
@with_kw mutable struct FirmCurrentAccount
    # Inflows
    S::Float64 = 0.0               # Sales
    add_debt::Float64 = 0.0        # Additional debts
    rev_dep::Float64 = 0.0         # Deposit revenues

    # Outflows
    TCL::Float64 = 0.0             # Total cost of labor
    TCI::Float64 = 0.0             # Total cost of investments
    TCE::Float64 = 0.0             # Total cost of energy
    int_debt::Float64 = 0.0        # Interest paid over debt
    rep_debt::Float64 = 0.0        # Debt repayments
end


"""
Clears firms current account for next period.
"""
function clear_firm_currentaccount_p!( 
    ca::FirmCurrentAccount
    )::FirmCurrentAccount

    ca.S = 0.0
    ca.add_debt = 0.0             
    ca.rev_dep = 0.0

    ca.TCL = 0.0
    ca.TCI = 0.0
    ca.TCE = 0.0
    ca.int_debt = 0.0
    ca.rep_debt = 0.0

    return ca
end