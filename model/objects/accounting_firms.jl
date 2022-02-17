"""
BALANCE SHEET
"""
mutable struct Balance
    # Assets
    N :: Float64            # Inventories
    K :: Float64            # Capital
    NW :: Float64           # Liquid assets

    # Liabilities
    Deb :: Float64          # Debt
    EQ :: Float64           # Equity
end


"""
Closes balance by computing firm's equity.
If firm is insolvent, liquidate firm.
"""
function close_balance_p!(
    p::AbstractAgent,
    Λ_max::Float64
    )

    # TODO: make sure credit is paid off before investments!

    # Check if NW is negative, if so, get short-term credit
    if p.balance.NW < 0
        poss_credit = p.D[end] * Λ_max - p.balance.Deb
        # Check if difference can be borrowed
        if p.balance.NW < poss_credit
            # Firm borrows credit to finance shortage
            req_credit = -p.balance.NW
            borrow_funds_p!(p, req_credit)
        else
            # Firm can no longer finance expenses and should be liquidated
            # TODO
        end
    end

    if p == ConsumerGoodProducer
        p.balance.N = p.N_goods * p.p[end]
    end

    tot_assets = p.p[end] * p.balance.N + p.balance.K + p.balance.NW
    equity = tot_assets - p.balance.Deb
    p.balance.EQ = equity
end


"""
CURRENT ACCOUNT
"""
mutable struct FirmCurrentAccount
    # Inflows
    S :: Float64            # Sales
    Rev_Dep :: Float64      # Deposit revenues

    # Outflows
    TCL::Float64            # Total cost of labor
    TCI::Float64            # Total cost of investments
    int_Deb::Float64        # Interest paid over debt
    rep_Deb::Float64        # Debt repayments
end


"""
Clears firms current account for next period.
"""
function clear_firm_currentaccount!(
    ca::FirmCurrentAccount
    )

    ca.S = 0
    ca.Rev_Dep = 0
    ca.TCL = 0
    ca.TCI = 0
    ca.int_Deb = 0
    ca.rep_Deb = 0
end