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
    globalparam::GlobalParam,
    government,
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
        update_interest_payment_p!(model[p_id], globalparam.r)

        # Repay debts of period
        monthly_debt_payback_p!(model[p_id], globalparam.b)

        # Update valuation of capital stock
        writeoffs = update_K_p!(model[p_id], globalparam.η)

        # Compute profits
        compute_Π_p!(model[p_id], government; writeoffs)

        if t >= 1
            # Compute and pay energy and carbon taxes
            pay_carbontax_p!(model[p_id], government, t)
            pay_energytax_p!(model[p_id], government, t)
        end

        # Update liquid assets NW
        update_NW_p!(model[p_id], government.τᴾ)

        # Update valuation of inventory. kp do not have inventory.
        if typeof(model[p_id]) == ConsumerGoodProducer
            model[p_id].balance.N = model[p_id].p[end] * model[p_id].N_goods
        end

        # Determine how much the firm can have as NW at most
        max_NW = globalparam.max_NW_ratio * (model[p_id].curracc.TCL + model[p_id].curracc.TCI + 
                      model[p_id].curracc.int_debt + model[p_id].debt_installments[2])

        if typeof(model[p_id]) == ConsumerGoodProducer
            all_max_NW += max_NW
        end

        
        if model[p_id].balance.NW < 0

            # If not enough liquid assets available, borrow additional funds.
            borrow_funds_p!(model[p_id], -model[p_id].balance.NW, globalparam.b)
            model[p_id].balance.NW = 0

        elseif (!check_if_bankrupt_p!(model[p_id],  globalparam.t_wait) 
                && (model[p_id].balance.NW > max_NW) && (t > globalparam.t_wait))
            
            # If enough liquid assets available, first pay off debts and then pay out dividends
            excess_NW = model[p_id].balance.NW - max_NW

            # First pay off debts with excess funds, excess dividends are paid out
            firm_dividends = singular_debt_payback_p!(model[p_id], excess_NW)
            total_dividends += firm_dividends

            model[p_id].balance.NW = max_NW
            
        end

        # Compute Equity
        tot_assets = model[p_id].balance.N + model[p_id].balance.K + model[p_id].balance.NW
        model[p_id].balance.EQ = tot_assets - model[p_id].balance.debt

        if typeof(model[p_id]) == ConsumerGoodProducer
            all_NW += model[p_id].balance.NW
        end
    end

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
    profittax = p.curracc.profittax
    carbontax = p.curracc.carbontax
    energytax = p.curracc.energytax

    new_NW = (NW + S + add_debt + rev_dep - TCL - TCI - TCE - rep_debt - int_debt 
              - profittax - carbontax - energytax)

    if typeof(p) == ConsumerGoodProducer
        p.NW_growth = (new_NW - p.balance.NW) / p.balance.NW
    end
    
    p.balance.NW = new_NW

    p.true_c = p.Q[end] > 0 ? (TCL + TCI + TCE + rep_debt + add_debt + int_debt) / p.Q[end] : p.c[end]
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
Computes profit of p in previous time period and profit tax that has to be paid.
"""
function compute_Π_p!(
    p::AbstractAgent,
    government;
    writeoffs=0.0::Float64
    )
 
    Π = (p.curracc.S + p.curracc.rev_dep - p.curracc.TCL - p.curracc.TCE 
         - p.curracc.carbontax - p.curracc.energytax - p.curracc.int_debt - writeoffs)
    shift_and_append!(p.Π, Π)

    # Compute profit taxes
    p.curracc.profittax = max(government.τᴾ * p.Π[end], 0)
end


"""
Computes carbon tax p has to pay
"""
function pay_carbontax_p!(
    p::AbstractAgent,
    government,
    t::Int64
    )

    p.curracc.carbontax = government.τᶜ * p.emissions
    receive_carbontax_gov!(government, p.curracc.carbontax, t)
end


"""
Computes carbon tax p has to pay
"""
function pay_energytax_p!(
    p::AbstractAgent,
    government,
    t::Int64
    )

    p.curracc.energytax = government.τᴱ * p.EU
    receive_energytax_gov!(government, p.curracc.energytax, t)
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

    profittax::Float64 = 0.0       # Profit tax
    carbontax::Float64 = 0.0       # Carbon tax
    energytax::Float64 = 0.0       # Energy tax
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

    ca.profittax = 0.0
    ca.carbontax = 0.0
    ca.energytax = 0.0

    return ca
end