@Base.kwdef mutable struct Government <: AbstractAgent
    T::Int = T
    UB::Float64 = 70.0                      # unemployment benefits
    τᴵ::Float64 = 0.35                      # income tax
    τˢ::Float64 = 0.0                       # sales tax
    τᴾ::Float64 = 0.0                      # profit tax
    τᴱ::Float64 = 0.0                       # energy tax
    τᶜ::Float64 = 0.0                       # emission tax
    MS::Float64 = 0.0                       # money stock owned by government
    curracc::GovCurrentAccount = GovCurrentAccount(T=T) # current account of government spending
end


# function initialize_government()
#     gov_struct = Government(
#         70,                          # UB: unemployment benefits
#         0.35,                           # τᴵ: income tax
#         0.0,                            # τˢ: sales tax
#         0.35,                           # τᴾ: profit tax
#         0.0,                            # τᴱ: energy tax
#         0.0,                            # τᶜ: emission tax
#         0.0,                            # MS: money stock owned by government
#         GovCurrentAccount(              # curr_account: current account of government spending
#             [],
#             [],
#             [],
#             [],
#             [],
#             [],
#             []
#         )  
#     )
#     return gov_struct
# end


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
        total_τᴵ += model[hh_id].I[end] * gov_struct.τᴵ
        push!(model[hh_id].Iᵀ, model[hh_id].I[end] * (1 - gov_struct.τᴵ))
    end

    # add total income tax to government current account
    # push!(gov_struct.curracc.Rev_τᴵ, total_τᴵ)
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
            p.Π[end] = p.Π[end] * (1 - gov_struct.τᴾ)
        end
    end

    # push!(gov_struct.curracc.Rev_τᴾ, total_τᴾ)
    gov_struct.curracc.Rev_τᴾ[t] = total_τᴾ
end


function compute_budget_balance(
    gov_struct::Government,
    t::Int
    )

    # TODO add sales, energy and emission tax
    # Tot_rev = ca.Rev_τᴵ[end] + ca.Rev_τˢ[end] + ca.Rev_τᴾ[end]
    Tot_rev = gov_struct.curracc.Rev_τᴵ[t] + gov_struct.curracc.Rev_τᴾ[t]

    # TODO add subsidies
    Tot_exp = gov_struct.curracc.Exp_UB[t]

    gov_balance = Tot_rev - Tot_exp
    # println("Gov deficit: ", gov_balance)

    # Pay off part of debt in case of positive balance
    gov_struct.MS += gov_balance
end


# function set_salestax_zero_gov!(
#     gov_struct::Government
#     )

#     push!(gov_struct.curracc.Rev_τˢ, 0)
# end

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

        total_I = sum(map(hh_id -> model[hh_id].I[end], all_hh))
        for hh_id in all_hh
            model[hh_id].W[end] += (model[hh_id].I[end] / total_I) * gov_struct.MS
        end
        gov_struct.MS = 0.0
    end
end