@Base.kwdef mutable struct Government <: AbstractAgent
    
    w_min::Float64 = 0.5
    UB::Float64 = 100 * 0.7 * w_min        # unemployment benefits

    # Tax rates
    τᴵ::Float64 = 0.25                     # income tax
    τᴷ::Float64 = 0.25                     # capital gains tax
    τˢ::Float64 = 0.                       # sales tax
    τᴾ::Float64 = 0.25                     # profit tax
    τᴱ::Float64 = 0.1                      # energy tax
    τᶜ::Float64 = 0.                       # emission tax

    MS::Float64 = 0.0                      # money stock owned by government
    curracc::GovCurrentAccount             # current account of government spending
    changedtaxrates::Union{Vector, Nothing} # vector of tax rates that will be changed once warmup period ends
end


function initgovernment(
    T::Int64,
    changedtaxrates::Union{Vector, Nothing}
    )::Government

    # Initialize government structure
    government = Government(curracc=GovCurrentAccount(T=T), changedtaxrates=changedtaxrates)

    return government
end


"""
Instates changed taxes at the end of the warmup period
"""
function instatetaxes!(
    government::Government
    )

    # If changed tax rates passed, change in government struct
    if government.changedtaxrates ≠ nothing
        for (taxtype, taxrate) in government.changedtaxrates
            setproperty!(government, taxtype, taxrate)
        end
    end
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
        receiveincome_hh!(model[hh_id], government.UB; isUB=true)
        total_UB += government.UB
    end

    # Add total UB spending to government current account
    government.curracc.Exp_UB[t] = total_UB
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
            shift_and_append!(p.Πᵀ, p.Π[end] * (1 - government.τᴾ))
        else
            shift_and_append!(p.Πᵀ, 0.0)
        end
    end

    government.curracc.Rev_τᴾ[t] = total_τᴾ
end


"""
Lets government receive income tax from employers
"""
function receive_incometax_gov!(
    government::Government,
    incometax::Float64,
    t::Int
    )

    government.curracc.Rev_τᴵ[t] += incometax
end


"""
Lets government receive capital gains tax from the indexfund
"""
function receive_capgains_tax_gov!(
    government::Government, 
    capgainstax::Float64,
    t::Int
    )

    government.curracc.Rev_τᴷ[t] += capgainstax
end


"""
Lets government receive sales taxes from consumer good producers
"""
function receive_salestax_gov!(
    government::Government,
    salestax::Float64,
    t::Int
    )

    government.curracc.Rev_τˢ[t] += salestax
end


"""
Lets government receive energy taxes from producers
"""
function receive_energytax_gov!(
    government::Government,
    energytax::Float64,
    t::Int
    )

    government.curracc.Rev_τᴱ[t] += energytax
end


"""
Lets government receive carbon taxes from producers
"""
function receive_carbontax_gov!(
    government::Government,
    carbontax::Float64,
    t::Int
    )

    government.curracc.Rev_τᶜ[t] += carbontax
end


function compute_budget_balance(
    government::Government,
    t::Int
    )

    Rev_τᴷ = t > 1 ? government.curracc.Rev_τᴷ[t-1] : 0.0

    # Compute total tax revenues
    Tot_rev = (government.curracc.Rev_τᴵ[t] + government.curracc.Rev_τᴾ[t] 
               + government.curracc.Rev_τˢ[t] + government.curracc.Rev_τᴱ[t] 
               + government.curracc.Rev_τᶜ[t] + Rev_τᴷ)

    # Compute total expediture
    Tot_exp = government.curracc.Exp_UB[t] + government.curracc.Exp_Sub[t]

    # Pay off part of debt in case of positive balance
    government.MS += (Tot_rev - Tot_exp)
end


"""
    resolve_gov_balance()

Redistributes surplusses such that government does not acquire assets, or acquires
    capital such that the government does not incur debts on the long run (long-run
    budget deficits thus drain the capital market)
"""
function resolve_gov_balance!(
    government::Government,
    indexfund,
    globalparam::GlobalParam,
    all_hh::Vector{Int64},
    model::ABM
    )

    # println("gov budget balance is $(government.MS)")

    if government.MS >= 0.0

        total_I = sum(hh_id -> model[hh_id].total_I > 0 ? model[hh_id].total_I ^ -globalparam.prog : 0., all_hh)

        for hh_id in all_hh
            share = model[hh_id].total_I > 0 ? (model[hh_id].total_I ^ -globalparam.prog) / total_I : 0.
            socialbenefits = government.MS * share
            receiveincome_hh!(model[hh_id], socialbenefits; socben=true)
        end
        government.MS = 0.0
    else
        issuegovbonds(indexfund, -government.MS)
        government.MS = 0.0
    end
end


"""
Determines the rate of income tax that should be paid over the agent's income.
    This will be a flat rate by default, and can be ammended to be a progressive tax.
"""
function determine_incometaxrate(
    government::Government,
    income::Float64
    )::Float64

    return government.τᴵ
end