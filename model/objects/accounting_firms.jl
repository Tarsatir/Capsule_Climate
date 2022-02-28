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
function close_balance_cp!(
    cp::AbstractAgent,
    Λ::Float64,
    ΔD::Float64,
    r::Float64
    )

    # Repay debts of period
    payback_debt_p!(cp)

    # Update valuation of inventory
    cp.balance.N = cp.p[end] * cp.N_goods

    # Update valuation of capital stock
    cp.balance.K = sum(map(machine -> machine.freq * machine.p, cp.Ξ))
    # TODO: include write-offs

    # Compute balance of current account
    curracc_in = cp.curracc.S + cp.curracc.rev_dep
    curracc_out = cp.curracc.TCL + cp.curracc.TCI + cp.curracc.int_Deb + cp.curracc.rep_Deb
    curracc_balance = curracc_in - curracc_out
    # println(curracc_balance)

    # Decide if extra credit is needed and how much
    if curracc_balance < 0.0
        req_borrowing = 0
        if -curracc_balance > cp.balance.NW
            # Extra credit is needed, borrow extra funds
            curr_Deb = sum(cp.Deb_installments)
            max_Deb = Λ * cp.curracc.S
            req_borrowing = -curracc_balance - cp.balance.NW
            add_borrowing = max(min(max_Deb - curr_Deb, req_borrowing), 0)

            if add_borrowing > 0
                borrow_funds_p!(cp, add_borrowing)
            end
        else
            # Current account deficit can be financed from liquid assets
            # cp.balance.NW -= -curracc_balance
        end
    end
    
    # Borrow liquid assets if needed
    # if curracc_balance < 0
    #     borrow_funds_p!(cp, -curracc_balance)
    # end

    # Compute profits
    compute_Π_cp!(cp, r)

    # println(cp.Π[end], " ", cp.cI)

    cp.balance.Deb = sum(cp.Deb_installments)
    cp.balance.NW = cp.balance.NW + cp.Π[end] - cp.cI

    # Compute Equity
    tot_assets = cp.balance.N + cp.balance.K + cp.balance.NW
    cp.balance.EQ = tot_assets - cp.balance.Deb
    # println(tot_assets, " ", cp.balance.EQ)

    if cp.balance.EQ < 0
        # Liquidate firm
    end
end


"""
Closes balance by computing firm's equity.
If firm is insolvent, liquidate firm.
"""
function close_balance_kp!(
    kp::AbstractAgent,
    Λ::Float64,
    ΔD::Float64
    )

    kp.balance.Deb = sum(kp.Deb_installments)
    kp.balance.NW = kp.balance.NW + kp.Π[end] - kp.Deb_installments[1]

    # Compute Equity
    tot_assets = kp.balance.N + kp.balance.K + kp.balance.NW
    equity = tot_assets - kp.balance.Deb
    kp.balance.EQ = equity

    if kp.balance.EQ < 0
        # Liquidate firm
    end
end


"""
CURRENT ACCOUNT
"""
mutable struct FirmCurrentAccount
    # Inflows
    S :: Float64            # Sales
    rev_dep :: Float64      # Deposit revenues

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
    ca.rev_dep = 0
    ca.TCL = 0
    ca.TCI = 0
    ca.int_Deb = 0
    ca.rep_Deb = 0
end