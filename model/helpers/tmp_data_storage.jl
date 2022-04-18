"""
Mutable struct that holds the data structures used to save data from the consumer
    market process.
"""
@with_kw mutable struct CMData
    n_hh::Int
    n_cp::Int
    all_C::Vector{Float64} = zeros(Float64, n_hh)
    all_N::Vector{Float64} = zeros(Float64, n_cp)
    sold_per_hh::Vector{Float64} = zeros(Float64, n_hh)
    sold_per_cp::Vector{Float64} = zeros(Float64, n_cp)
    weights::Matrix{Float64} = rand(n_hh, n_cp)
    weights_sum::Vector{Float64} = zeros(Float64, n_hh)
    transactions::Matrix{Float64} = zeros(Float64, n_hh, n_cp)
    frac_sellable::Vector{Float64} = ones(Float64, n_cp)
    C_spread::Matrix{Float64} = zeros(Float64, n_hh, n_cp)
    demand_per_cp::Vector{Float64} = zeros(Float64, n_cp)
end