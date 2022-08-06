Base.@kwdef mutable struct IndexFund
    T::Int=T 
    Assets::Float64 = 0.0                           # Total amount of assets alocated in the fund
    funds_inv::Vector{Float64} = zeros(Float64, T)  # Amount of funds used for investment in last period
end


"""
Initializes index fund struct
"""
function initialize_indexfund()::IndexFund

    indexfund = IndexFund()
    
    return indexfund
end


"""
Takes dividends from producers
"""
function receive_dividends_if!(
    indexfund::IndexFund,
    dividends::Float64
    )
    
    indexfund.Assets += dividends
end


"""
Deducts funds for investment in company
"""
function decide_investments_if!(
    indexfund::IndexFund,
    all_req_NW::Float64,
    t::Int
    )::Float64

    frac_NW_if = indexfund.Assets > 0 ? min(indexfund.Assets / all_req_NW, 1.0) : 0.0
    indexfund.funds_inv[t] = all_req_NW * frac_NW_if
    indexfund.Assets -= indexfund.funds_inv[t]
    return frac_NW_if
end


"""
Deducts funds for net debts lost
"""
function deduct_unpaid_net_debts_if!(
    indexfund::IndexFund,
    total_unpaid_net_debt::Float64
    )

    indexfund.Assets -= total_unpaid_net_debt
end


"""
Distributes dividends over participants in indexfund
"""
function distribute_dividends_if!(
    indexfund::IndexFund,
    government::Government,
    all_hh::Vector{Int64},
    τᴷ::Float64,
    t::Int64,
    model::ABM
    )

    # Distribute proportional to wealth
    all_W = map(hh_id -> model[hh_id].W, all_hh)
    total_wealth = sum(all_W)

    total_capgains_tax = 0

    for hh_id in all_hh
        # Do not award to most wealthy household to avoid one household
        # taking over all wealth
        if t < 10 
            dividend_share = (model[hh_id].W / total_wealth) * indexfund.Assets
        elseif model[hh_id].W ≠ maximum(all_W)
            dividend_share = (model[hh_id].W / (total_wealth - maximum(all_W))) * indexfund.Assets
        else
            dividend_share = 0.
        end

        total_capgains_tax += τᴷ * dividend_share
        receiveincome_hh!(model[hh_id], dividend_share * (1 - τᴷ); capgains=true)
    end

    # Reset assets back to zero
    indexfund.Assets = 0.

    receive_capgains_tax_gov!(government, total_capgains_tax, t)
end


"""
Lets government issue bonds on the capital market.
"""
function issuegovbonds(
    indexfund::IndexFund, 
    govdeficit::Float64
    )

    indexfund.Assets -= govdeficit
end