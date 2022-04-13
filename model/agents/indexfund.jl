Base.@kwdef mutable struct IndexFund
    T::Int=T 
    Assets::Float64 = 0.0                           # Total amount of assets alocated in the fund
    funds_inv::Vector{Float64} = zeros(Float64, T)  # Amount of funds used for investment in last period
end


"""
Initializes index fund struct
"""
function initialize_indexfund()::IndexFund

    indexfund_struct = IndexFund()
    
    return indexfund_struct
end


"""
Takes dividends from producers
"""
function receive_dividends_if!(
    indexfund_struct::IndexFund,
    dividends::Float64
    )
    
    indexfund_struct.Assets += dividends
end


"""
Deducts funds for investment in company
"""
function decide_investments_if!(
    indexfund_struct::IndexFund,
    all_req_NW::Float64,
    t::Int
    )::Float64

    frac_NW_if = max(0, min(indexfund_struct.Assets / all_req_NW, 1))
    indexfund_struct.funds_inv[t] = all_req_NW * frac_NW_if
    indexfund_struct.Assets -= indexfund_struct.funds_inv[t]
    return frac_NW_if
end


"""
Deducts funds for net debts lost
"""
function deduct_unpaid_net_debts_if!(
    indexfund_struct::IndexFund,
    total_unpaid_net_debt::Float64
    )
    indexfund_struct.Assets -= total_unpaid_net_debt
end


"""
Distributes dividends over participants in indexfund
"""
function distribute_dividends_if!(
    indexfund_struct::IndexFund,
    all_hh::Vector{Int},
    model::ABM
    )

    # Distribute proportional to wealth
    total_wealth = sum(hh_id -> model[hh_id].W[end], all_hh)

    for hh_id in all_hh

        dividend_share = (model[hh_id].W[end] / total_wealth) * indexfund_struct.Assets

        # TODO: make separate function in hh for this
        model[hh_id].W[end] += dividend_share
    end

    # Reset assets back to zero
    indexfund_struct.Assets = 0.0

end