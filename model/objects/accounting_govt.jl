mutable struct GovCurrentAccount
    # Revenues
    Rev_τᴵ :: Vector{Float32}            # hist revenues of income tax
    Rev_τˢ :: Vector{Float32}            # hist revenues of sales tax
    Rev_τᴾ :: Vector{Float32}            # hist revenues of profit tax
    Rev_τᴱ :: Vector{Float32}            # hist revenues of energy tax
    Rev_τᶜ :: Vector{Float32}            # hist revenues of emission tax

    # Expenditures
    Exp_UB :: Vector{Float32}            # hist expenditures on unemployment benefits
    Exp_Sub :: Vector{Float32}           # hist expenditures on subsidies

end