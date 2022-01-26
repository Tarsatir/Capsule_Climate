
mutable struct Government <: AbstractAgent
    UB :: Float64                       # unemployment benefits
    τᴵ :: Float64                       # income tax
    τˢ :: Float64                       # sales tax
    τᴾ :: Float64                       # profit tax
    τᴱ :: Float64                       # energy tax
    τᶜ :: Float64                       # emission tax
end

function initialize_government()
    gov_struct = Government(
        50,                             # UB: unemployment benefits
        0.2,                            # τᴵ: income tax
        0.2,                            # τˢ: sales tax
        0.2,                            # τᴾ: profit tax
        0.0,                            # τᴱ: energy tax
        0.0                             # τᶜ: emission tax
    )
    return gov_struct
end

function pay_unemployment_benefits_gov!(gov_struct, unemployed)

    for hh in unemployed
        push!(hh.I, gov_struct.UB)
    end

end