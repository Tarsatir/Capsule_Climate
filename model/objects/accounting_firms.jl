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
    Λ::Float64,
    ΔD::Float64
    )

    # println(p.curracc, " ", p.Deb_installments)

    p.balance.Deb = sum(p.Deb_installments)

    # Check if additional credit is needed
    # if p.balance.NW < 0
    #     borrow_funds_p!(p, -p.balance.NW)
    # end

    # Repay debt (as far as NW allows)
    # Deb_to_be_repaid = min(p.balance.Deb * ΔD, p.balance.NW)
    # payback_debt_p!(p, Deb_to_be_repaid)

    # Compute Equity
    tot_assets = p.p[end] * p.balance.N + p.balance.K + p.balance.NW
    equity = tot_assets - p.balance.Deb
    p.balance.EQ = equity

    if p.balance.EQ < 0
        # Liquidate firm
    end
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
function clear_firm_currentaccount_p!(
    ca::FirmCurrentAccount
    )

    ca.S = 0
    ca.Rev_Dep = 0
    ca.TCL = 0
    ca.TCI = 0
    ca.int_Deb = 0
    ca.rep_Deb = 0
end