"""
BALANCE SHEET
"""
mutable struct Balance
    # Assets
    N :: Float64            # Inventories
    K :: Float64            # Capital
    NW :: Float64           # Liquid assets

    # Liabilities
    debt :: Float64         # Debt
    EQ :: Float64           # Equity
end


"""
Closes balance by computing firm's equity.
If firm is insolvent, liquidate firm.
"""
function close_balance_p!(
    p::AbstractAgent,
    Λ::Float64,
    r::Float64,
    η::Int,
    τᴾ::Float64,
    indexfund_struct
    )

    p.balance.debt = max(0.0, p.balance.debt)

    # Compute interest payment
    update_interest_payment_p!(p, r)

    # Repay debts of period
    payback_debt_p!(p)

    # Update valuation of capital stock
    writeoffs = update_K_p!(p, η)

    # Compute profits
    compute_Π_p!(p, writeoffs)

    # Update liquid assets NW
    update_NW_p!(p, τᴾ)

    # If not enough liquid assets available, borrow additional funds.
    if p.balance.NW < 0
        max_add_debt = p.curracc.S * Λ - p.balance.debt
        add_debt = min(-p.balance.NW, max_add_debt)
        borrow_funds_p!(p, add_debt)
        p.balance.NW += add_debt
        p.balance.debt += add_debt
    end

    # Update valuation of inventory. kp do not have inventory.
    if typeof(p) == ConsumerGoodProducer
        p.balance.N = p.p[end] * p.N_goods
    end

    # Compute Equity
    tot_assets = p.balance.N + p.balance.K + p.balance.NW
    p.balance.EQ = tot_assets - p.balance.debt

    # TRIAL OF INDEX FUND
    # d = 2

    # max_net_NW = d * (p.curracc.TCL + p.curracc.TCI + p.curracc.int_debt + p.curracc.rep_debt)

    # if p.balance.NW - p.balance.debt > max_net_NW && !check_if_bankrupt_p!(p)
    #     indexfund_struct.Assets += (p.balance.NW - p.balance.debt - max_net_NW)
    #     p.balance.NW = max_net_NW + p.balance.debt
    #     tot_assets = p.balance.N + p.balance.K + p.balance.NW
    #     p.balance.EQ = tot_assets - p.balance.debt
    # end

    # If NW is negative, maximum debt is reached, and EQ is set to
    # a negative value so the firm is declared bankrupt
    if p.balance.NW < 0
        p.balance.EQ = -1.0
    end
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
        K = 0
        writeoffs = 0
        for machine in p.Ξ
            newval = machine.p * machine.freq
            K += max((η - machine.age) / η, 0) * newval
            if machine.age <= η
                writeoffs += (1 / η) * newval
            end
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
    rep_debt = p.curracc.rep_debt
    add_debt = p.curracc.add_debt
    int_debt = p.curracc.int_debt
    rev_dep = p.curracc.rev_dep
    profit_tax = max(τᴾ * p.Π[end], 0)

    # if typeof(p) == ConsumerGoodProducer
    #     println("NW: $NW, S: $S, TCL: $TCL, TCI: $TCI, rd: $rep_debt, ad: $add_debt, id: $int_debt, pt: $profit_tax")
    # end

    new_NW = NW + S + add_debt + rev_dep - TCL - TCI - rep_debt - int_debt - profit_tax

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
    p::AbstractAgent,
    writeoffs=0.0::Float64
    )
 
    Π = p.curracc.S + p.curracc.rev_dep - p.curracc.TCL - p.curracc.int_debt
    # if typeof(p) == CapitalGoodProducer
    #     println("profit ", Π, " S ", p.curracc.S, "  TCL ", p.curracc.TCL, " TCI ", p.curracc.TCI, " int ", p.curracc.int_debt, " p ", p.p[end])
    #     println("labor ", p.ΔLᵈ, " O ", p.O, " O/B ", p.O/p.B[end], " RD ", p.RD[end] / p.w̄[end], " L ", p.L, " w ", p.w̄[end])
    # end
    # Π = p.curracc.S + p.curracc.rev_dep - p.curracc.TCL - p.curracc.int_debt - writeoffs

    push!(p.Π, Π)
end


"""
CURRENT ACCOUNT
"""
mutable struct FirmCurrentAccount
    # Inflows
    S :: Float64            # Sales
    add_debt ::Float64       # Additional debts
    rev_dep :: Float64      # Deposit revenues

    # Outflows
    TCL::Float64            # Total cost of labor
    TCI::Float64            # Total cost of investments
    int_debt::Float64        # Interest paid over debt
    rep_debt::Float64        # Debt repayments
end


"""
Clears firms current account for next period.
"""
function clear_firm_currentaccount_p!(
    ca::FirmCurrentAccount
    )

    ca.S = 0
    ca.add_debt = 0
    ca.rev_dep = 0
    ca.TCL = 0
    ca.TCI = 0
    ca.int_debt = 0
    ca.rep_debt = 0
end