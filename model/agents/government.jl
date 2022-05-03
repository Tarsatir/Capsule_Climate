@Base.kwdef mutable struct Government <: AbstractAgent
    UB::Float64 = 0.0                     # unemployment benefits
    w_min::Float64 = 0.0
    # Tax rates
    τᴵ::Float64 = 0.0                      # income tax
    τˢ::Float64 = 0.0                      # sales tax
    τᴾ::Float64 = 0.0                      # profit tax
    τᴱ::Float64 = 0.0                      # energy tax
    τᶜ::Float64 = 0.0                      # emission tax

    MS::Float64 = 0.0                      # money stock owned by government
    curracc::GovCurrentAccount             # current account of government spending
end


"""
Loops over all unemployed households and pays UB
"""
function pay_unemployment_benefits_gov!(
    gov_struct, 
    unemployed::Vector{Int},
    t::Int,
    model::ABM
    )
    # pay out unemployment benefits to households
    total_UB = 0
    for hh_id in unemployed
        # push!(hh.I, gov_struct.UB)
        get_income_hh!(model[hh_id], gov_struct.UB)
        total_UB += gov_struct.UB
    end

    # add total UB spending to government current account
    # push!(gov_struct.curracc.Exp_UB, total_UB)
    gov_struct.curracc.Exp_UB[t] = total_UB
end


"""
Levies income tax on all households
"""
function levy_income_tax_gov!(
    gov_struct::Government, 
    all_hh::Vector{Int},
    t::Int,
    model::ABM
    )

    total_τᴵ = 0
    for hh_id in all_hh
        total_τᴵ += model[hh_id].I * gov_struct.τᴵ
        model[hh_id].Iᵀ = model[hh_id].I * (1 - gov_struct.τᴵ)
    end

    # add total income tax to government current account
    gov_struct.curracc.Rev_τᴵ[t] = total_τᴵ
end


"""
Levies profit tax on all producers
"""
function levy_profit_tax_gov!(
    gov_struct::Government,
    all_p::Vector{Int},
    t::Int,
    model::ABM
    )

    total_τᴾ = 0
    for p_id in all_p
        p = model[p_id]
        # Only levy tax when profit is positive
        if p.Π[end] > 0
            total_τᴾ += p.Π[end] * gov_struct.τᴾ
            # p.Π[end] = p.Π[end] * (1 - gov_struct.τᴾ)
            shift_and_append!(p.Πᵀ, (1 - gov_struct.τᴾ))
        else
            shift_and_append!(p.Πᵀ, 0.0)
        end
    end

    gov_struct.curracc.Rev_τᴾ[t] = total_τᴾ
end


function compute_budget_balance(
    gov_struct::Government,
    t::Int
    )

    # Compute total tax revenues
    Tot_rev = (gov_struct.curracc.Rev_τᴵ[t] + gov_struct.curracc.Rev_τᴾ[t] + gov_struct.curracc.Rev_τˢ[t]
               + gov_struct.curracc.Rev_τᴱ[t] + gov_struct.curracc.Rev_τᶜ[t])

    # Compute total expediture
    Tot_exp = gov_struct.curracc.Exp_UB[t] + gov_struct.curracc.Exp_Sub[t]

    # Pay off part of debt in case of positive balance
    gov_struct.MS += (Tot_rev - Tot_exp)
end


function add_salestax_transaction_gov!(
    gov_struct::Government,
    sales_tax_bp::Float64, 
    sales_tax_lp::Float64,
    t::Int
    )
    gov_struct.curracc.Rev_τˢ[t] += (sales_tax_bp + sales_tax_lp)
end


"""
Redistributes surplusses such that government does not acquire assets.
"""
function redistribute_surplus_gov!(
    gov_struct::Government,
    all_hh::Vector{Int},
    model
    )

    # NOTE: now does not apply redistribution between income groups!
    # NOTE: redist goes both ways!

    if gov_struct.MS > 0.0

        total_I = sum(hh_id -> model[hh_id].I, all_hh)
        for hh_id in all_hh
            model[hh_id].W += (model[hh_id].I / total_I) * gov_struct.MS
            # model[hh_id].W += gov_struct.MS / length(all_hh)
        end
        gov_struct.MS = 0.0
    end
end