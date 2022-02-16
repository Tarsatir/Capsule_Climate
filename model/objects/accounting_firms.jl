mutable struct Balance
    # Assets
    N :: Float64            # inventories
    K :: Float64            # capital
    NW :: Float64           # liquid assets

    # Liabilities
    Deb :: Float64          # debt
    EQ :: Float64           # equity
end


"""
Closes balance by computing firm's equity.
"""
function close_balance!(
    p::AbstractAgent
    )

    tot_assets = p.p[end] * p.balance.N + p.balance.K + p.balance.NW
    equity = tot_assets - p.balance.Deb
    p.balance.EQ = equity
end