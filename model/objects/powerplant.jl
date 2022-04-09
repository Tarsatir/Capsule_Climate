@with_kw mutable struct PowerPlant
    type::String                                # Type of the power plant ("Green" or "Dirty")
    c::Float64                                  # Marginal cost of production
    freq::Float64                               # Absolute frequency of machine
    capacity::Float64                           # Capacity to produce energy
    Aᵗ::Float64                                 # Thermal efficiency
    em::Float64                                 # emissions per unit of energy
end


"""
Updates the marginal production cost of the power plant based on the type and 
    price of fossil fuels
    Lamperti (2018), eq 11.
"""
function update_c_pp!(
    pp::PowerPlant, 
    p_f::Float64
    )

    pp.c = pp.type == "Dirty" ? p_f / pp.Aᵗ : 0.0
end