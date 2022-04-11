@with_kw mutable struct Climate

    T::Int
    
    carbon_emissions::Vector{Float64} = zeros(Float64, T)
end