@Base.kwdef mutable struct Government <: AbstractAgent
    UB::Float64 = 40.0                     # unemployment benefits
    w_min::Float64 = 0.5

    # Tax rates
    τᴵ::Float64 = 0.3                      # income tax
    τᴷ::Float64 = 0.05                     # capital gains tax
    τˢ::Float64 = 0.0                      # sales tax
    τᴾ::Float64 = 0.3                      # profit tax
    τᴱ::Float64 = 0.0                      # energy tax
    τᶜ::Float64 = 0.0                      # emission tax

    MS::Float64 = 0.0                      # money stock owned by government
    curracc::GovCurrentAccount             # current account of government spending
end


"""
Loops over all unemployed households and pays UB
"""
function pay_unemployment_benefits_gov!(
    government::Government, 
    unemployed::Vector{Int},
    t::Int,
    model::ABM
    )

    # Pay out unemployment benefits to households
    total_UB = 0
    for hh_id in unemployed
        receiveincome_hh!(model[hh_id], government.UB)
        total_UB += government.UB
    end

    # Add total UB spending to government current account
    government.curracc.Exp_UB[t] = total_UB
end


"""
Levies income tax on all households
"""
function levy_income_tax_gov!(
    government::Government, 
    all_hh::Vector{Int},
    t::Int,
    model::ABM
    )

    total_τᴵ = 0
    for hh_id in all_hh
        total_τᴵ += model[hh_id].I * government.τᴵ
        model[hh_id].Iᵀ = model[hh_id].I * (1 - government.τᴵ)
    end

    # Add total income tax to government current account
    government.curracc.Rev_τᴵ[t] = total_τᴵ
end


"""
Levies profit tax on all producers
"""
function levy_profit_tax_gov!(
    government::Government,
    all_p::Vector{Int},
    t::Int,
    model::ABM
    )

    total_τᴾ = 0
    for p_id in all_p
        p = model[p_id]
        # Only levy tax when profit is positive
        if p.Π[end] > 0
            total_τᴾ += p.Π[end] * government.τᴾ
            shift_and_append!(p.Πᵀ, (1 - government.τᴾ))
        else
            shift_and_append!(p.Πᵀ, 0.0)
        end
    end

    government.curracc.Rev_τᴾ[t] = total_τᴾ
end


"""
Lets government receive capital gains tax from the indexfund
"""
function receive_capgains_tax_gov!(
    government::Government, 
    total_capgains_tax::Float64,
    t::Int
    )

    government.curracc.Rev_τᴷ[t] += total_capgains_tax
    government.MS += total_capgains_tax
end


function compute_budget_balance(
    government::Government,
    t::Int
    )

    # Compute total tax revenues
    Tot_rev = (government.curracc.Rev_τᴵ[t] + government.curracc.Rev_τᴷ[t] 
               + government.curracc.Rev_τᴾ[t] + government.curracc.Rev_τˢ[t]
               + government.curracc.Rev_τᴱ[t] + government.curracc.Rev_τᶜ[t])

    # Compute total expediture
    Tot_exp = government.curracc.Exp_UB[t] + government.curracc.Exp_Sub[t]

    # Pay off part of debt in case of positive balance
    government.MS += (Tot_rev - Tot_exp)
end


function add_salestax_transaction_gov!(
    government::Government,
    sales_tax_bp::Float64, 
    sales_tax_lp::Float64,
    t::Int
    )
    government.curracc.Rev_τˢ[t] += (sales_tax_bp + sales_tax_lp)
end


"""
Redistributes surplusses such that government does not acquire assets.
"""
function redistribute_surplus_gov!(
    government::Government,
    all_hh::Vector{Int},
    model
    )

    # NOTE: now does not apply redistribution between income groups!
    # NOTE: redist goes both ways!

    if government.MS > 0.0

        total_I = sum(hh_id -> model[hh_id].I, all_hh)
        for hh_id in all_hh
            model[hh_id].W += (model[hh_id].I / total_I) * government.MS
            # model[hh_id].W += government.MS / length(all_hh)
        end
        government.MS = 0.0
    end
end