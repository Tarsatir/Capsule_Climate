mutable struct Bank (
    Deposits::Vector{Tuple{Int, Float64}}               # deposits at bank
    Debts::Vector{Tuple{Int, Float64}}                  # debts at bank
)