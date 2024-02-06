@with_kw mutable struct GovCurrentAccount

    T::Int=T

    # Revenues
    Rev_τᴵ::Vector{Float64} = zeros(Float64, T)  # hist revenues of income tax
    Rev_τᴷ::Vector{Float64} = zeros(Float64, T)  # hist revenues of capital gains tax
    Rev_τˢ::Vector{Float64} = zeros(Float64, T)  # hist revenues of sales tax
    Rev_τᴾ::Vector{Float64} = zeros(Float64, T)  # hist revenues of profit tax
    Rev_τᴱ::Vector{Float64} = zeros(Float64, T)  # hist revenues of energy tax
    Rev_τᶜ::Vector{Float64} = zeros(Float64, T)  # hist revenues of emission tax

    # Expenditures
    Exp_UB::Vector{Float64} = zeros(Float64, T)  # hist expenditures on unemployment benefits
    Exp_Sub::Vector{Float64} = zeros(Float64, T) # hist expenditures on subsidies
end