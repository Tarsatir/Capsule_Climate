mutable struct Government <: AbstractAgent
    UB :: Float64                       # unemployment benefits
    τᴵ :: Float64                       # income tax
    τˢ :: Float64                       # sales tax
    τᴾ :: Float64                       # profit tax
    τᴱ :: Float64                       # energy tax
    τᶜ :: Float64                       # emission tax
    curr_acc :: GovCurrentAccount       # current account of government spending
end


function initialize_government()
    gov_struct = Government(
        50,                             # UB: unemployment benefits
        0.2,                            # τᴵ: income tax
        0.2,                            # τˢ: sales tax
        0.2,                            # τᴾ: profit tax
        0.0,                            # τᴱ: energy tax
        0.0,                             # τᶜ: emission tax
        initialize_govcurrentaccount()  # curr_account: current account of government spending
    )
    return gov_struct
end


function pay_unemployment_benefits_gov!(gov_struct, unemployed)

    # pay out unemployment benefits to households
    total_UB = 0
    for hh in unemployed
        push!(hh.I, gov_struct.UB)
        total_UB += gov_struct.UB
    end

    # add total UB spending to government current account
    push!(gov_struct.curr_acc.Exp_UB, total_UB)
end


"""
Levies income tax on all households
"""
function levy_income_tax_gov!(gov_struct, all_hh)

    total_τᴵ = 0
    for hh in all_hh
        total_τᴵ += hh.I[end] * gov_struct.τᴵ
        hh.I[end] = hh.I[end] * (1 - gov_struct.τᴵ)
    end

    # add total income tax to government current account
    push!(gov_struct.curr_acc.Rev_τᴵ, total_τᴵ)
end


function compute_budget_balance(gov_struct)

    println(gov_struct.curr_acc.Rev_τᴵ)
    println(gov_struct.curr_acc.Exp_UB)

    # Rev_tot = (gov_struct.curr_acc.Rev_τᴵ[end] + gov_struct.curr_acc.Rev_τˢ[end] + 
    #            gov_struct.curr_acc.Rev_τᴾ[end] + gov_struct.curr_acc.Rev_τᴱ[end] +
    #            gov_struct.curr_acc.Rev_τᶜ[end])
    
    # Exp_tot = gov_struct.curr_acc.Exp_UB[end] + gov_struct.curr_acc.Exp_Sub[end]

    # balance = Rev_tot - Exp_tot

    # println("gov balance: ", balance)

end