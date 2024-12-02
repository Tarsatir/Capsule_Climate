@Base.kwdef mutable struct Government <: AbstractAgent

    T::Int64
    w_min::Float64 = 0.5
    UB::Float64 = 100 * 0.7 * w_min        # unemployment benefits

    # (Default) tax rates
    τᴵ::Float64 = 0.25                              # income tax
    τᴷ::Float64 = 0.25                              # capital gains tax
    τˢ::Float64 = 0.                                # sales tax
    τᴾ::Float64 = 0.25                              # profit tax
    τᴱ::Float64 = 0.                                # energy tax
    τᶜ::Float64 = 0.                                # emission tax

    # Tax rates over time
    τᴵ_ts::Vector{Float64} = @MVector fill(τᴵ, T)   # income tax
    τᴷ_ts::Vector{Float64} = @MVector fill(τᴷ, T)   # capital gains tax
    τˢ_ts::Vector{Float64} = @MVector fill(τˢ, T)   # sales tax
    τᴾ_ts::Vector{Float64} = @MVector fill(τᴾ, T)   # profit tax
    τᴱ_ts::Vector{Float64} = @MVector fill(τᴱ, T)   # energy tax
    τᶜ_ts::Vector{Float64} = @MVector fill(τᶜ, T)   # emission tax


    MS::Float64 = 0.0                      # money stock owned by government

    # Revenues
    rev_incometax::Vector{Float64} = zeros(Float64, T)  # hist revenues of income tax
    rev_capitaltax::Vector{Float64} = zeros(Float64, T)  # hist revenues of capital gains tax
    rev_salestax::Vector{Float64} = zeros(Float64, T)  # hist revenues of sales tax
    rev_profittax::Vector{Float64} = zeros(Float64, T)  # hist revenues of profit tax
    rev_energytax::Vector{Float64} = zeros(Float64, T)  # hist revenues of energy tax
    rev_carbontax::Vector{Float64} = zeros(Float64, T)  # hist revenues of emission tax

    # Expenditures
    exp_UB::Vector{Float64} = zeros(Float64, T)  # hist expenditures on unemployment benefits
    exp_subsidies::Vector{Float64} = zeros(Float64, T) # hist expenditures on subsidies

    # changed_taxrates::Union{Vector, Nothing} # vector of tax rates that will be changed once warmup period ends
end


function initgovernment(
    T::Int64,
    t_warmup::Int64,
    changed_taxrates::Union{Vector, Nothing}
)::Government

    government = Government(T=T)

    if !isnothing(changed_taxrates)
        for changed_taxrate in changed_taxrates
            change_taxrate!(changed_taxrate, government, t_warmup)
        end
    end

    return government
end


function change_taxrate!(
    changed_taxrate::Tuple, 
    government::Government,
    t_warmup::Int64
)
    taxtype = changed_taxrate[1]
    taxtype_ts = Symbol(String(taxtype) * "_ts")
    original_taxrate = getproperty(government, taxtype)

    # if length(changed_taxrate) == 2
    #     # Change to new tax rate at end of warmup period
    #     taxrate = changed_taxrate[2]
    #     tax_ts = getproperty(government, taxtype_ts)
    #     tax_ts[t_warmup:end] .= taxrate
    #     #in case of shock experiment
    #     # change return taxrate to inital value after 100 entries
    #     tax_ts[t_warmup+20:end] .= original_taxrate
    if length(changed_taxrate) == 2
        # Change to new tax rate at the end of warmup period
        taxrate = changed_taxrate[2]
        tax_ts = getproperty(government, taxtype_ts)
        tax_ts[t_warmup:end] .= taxrate
    else
        # Linear interpolation between τ_start and τ_end
        taxrate_start = changed_taxrate[2]
        taxrate_end = changed_taxrate[3]
        tax_ts = getproperty(government, Symbol(String(taxtype) * "_ts"))
        tax_ts[begin:t_warmup] .= taxrate_start 
        tax_ts[t_warmup:end] .= Base._linspace(taxrate_start, taxrate_end, length(tax_ts) - t_warmup + 1)
    end

    setproperty!(government, taxtype_ts, tax_ts)
end


# function change_taxrate!(
#     changed_taxrate::Tuple{Symbol, Float64, Float64}, 
#     government::Government,
#     t_warmup::Int64
# )
#     taxtype = changed_taxrate[1]
#     taxtype_ts = Symbol(String(taxtype) * "_ts")

#     taxrate_start = changed_taxrate[2]
#     taxrate_end = changed_taxrate[3]
#     tax_ts = getproperty(government, Symbol(String(taxtype) * "_ts"))
#     tax_ts[begin:t_warmup] .= taxrate_start
#     tax_ts[t_warmup:end] .= Base._linspace(taxrate_start, taxrate_end, length(tax_ts) - t_warmup + 1)

#     setproperty!(government, taxtype_ts, tax_ts)
# end


# """
# Instates changed taxes at the end of the warmup period
# """
# function instate_taxes!(
#     government::Government,
#     t::Int64,
#     t_warmup::Int64
# )

#     # If changed tax rates passed, change in government struct
#     if !isnothing(government.changed_taxrates)
#         for (taxtype, taxrate) in government.changed_taxrates
#             setproperty!(government, taxtype, taxrate)
#         end
#     end
# end


function update_taxrates!(
    government::Government,
    t::Int64
)
    government.τᴵ = government.τᴵ_ts[t]
    government.τᴷ = government.τᴷ_ts[t]
    government.τˢ = government.τˢ_ts[t]
    government.τᴾ = government.τᴾ_ts[t]
    government.τᴱ = government.τᴱ_ts[t]
    government.τᶜ = government.τᶜ_ts[t]
end


"""
Loops over all unemployed households and pays UB
"""
function pay_unemployment_benefits_gov!(
    government::Government, 
    unemployed::Vector{Int},
    t::Int64,
    model::ABM
)

    # Pay out unemployment benefits to households
    total_UB = 0
    for hh_id in unemployed
        receiveincome_hh!(model[hh_id], government.UB; isUB=true)
        total_UB += government.UB
    end

    # Add total UB spending to government current account
    government.exp_UB[t] = total_UB
end


"""
Levies profit tax on all producers
"""
function levy_profit_tax_gov!(
    government::Government,
    all_p::Vector{Int64},
    t::Int64,
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

    government.rev_profittax[t] = total_τᴾ
end


"""
Lets government receive income tax from employers
"""
function receive_incometax_gov!(
    government::Government,
    incometax::Float64,
    t::Int64
)

    government.rev_incometax[t] += incometax
end


"""
Lets government receive capital gains tax from the indexfund
"""
function receive_capgains_tax_gov!(
    government::Government, 
    capgainstax::Float64,
    t::Int64
)
    government.rev_capitaltax[t] += capgainstax
end


"""
Lets government receive sales taxes from consumer good producers
"""
function receive_salestax_gov!(
    government::Government,
    salestax::Float64,
    t::Int64
)
    government.rev_salestax[t] += salestax
end


"""
Lets government receive energy taxes from producers
"""
function receive_energytax_gov!(
    government::Government,
    energytax::Float64,
    t::Int64
)
    government.rev_energytax[t] += energytax
end


"""
Lets government receive carbon taxes from producers
"""
function receive_carbontax_gov!(
    government::Government,
    carbontax::Float64,
    t::Int64
)
    government.rev_carbontax[t] += carbontax
end


function compute_budget_balance(
    government::Government,
    t::Int64
)

    rev_capitaltax = t > 1 ? government.rev_capitaltax[t-1] : 0.0

    # Compute total tax revenues
    total_revenue = (government.rev_incometax[t] + government.rev_profittax[t] 
                    + government.rev_salestax[t] + government.rev_energytax[t] 
                    + government.rev_carbontax[t] + rev_capitaltax)

    # Compute total expediture
    total_expenditure = government.exp_UB[t] + government.exp_subsidies[t]

    # Pay off part of debt in case of positive balance
    government.MS += (total_revenue - total_expenditure)
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
    t::Int64,
    model::ABM
)

    if government.MS >= 0.0

        total_I = sum(hh_id -> model[hh_id].total_I > 0 ? model[hh_id].total_I ^ -globalparam.prog : 0., all_hh)

        for hh_id in all_hh
            share = model[hh_id].total_I > 0 ? (model[hh_id].total_I ^ -globalparam.prog) / total_I : 0.
            socialbenefits = government.MS * share
            receiveincome_hh!(model[hh_id], socialbenefits; socben=true)
        end
        government.exp_subsidies[t] = government.MS
    else
        issuegovbonds(indexfund, -government.MS)
        government.MS = 0.0
    end
end


# """
# Determines the rate of income tax that should be paid over the agent's income.
#     This will be a flat rate by default, and can be ammended to be a progressive tax.
# """
# function determine_incometaxrate(
#     government::Government,
#     income::Float64
# )::Float64

#     return government.τᴵ
# end