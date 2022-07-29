@with_kw mutable struct PowerPlant
    type::String                                # Type of the power plant ("Green" or "Dirty")
    age::Int = 0                                # Age of power plant
    c::Float64                                  # Marginal cost of production
    freq::Float64                               # Absolute frequency of machine
    capacity::Float64                           # Capacity to produce energy
    Aᵀ::Float64                                 # Thermal efficiency
    em::Float64                                 # emissions per unit of energy
end

function init_powerplant(
    type::String,
    age::Int64,
    c::Float64,
    Aᵀ::Float64,
    emᵀ::Float64,
    globalparam::GlobalParam
    )

    return PowerPlant(
        type = type,
        age = age,
        c = c,
        freq = globalparam.freq_per_powerplant,
        capacity = type=="Dirty" ? globalparam.freq_per_powerplant * Aᵀ : globalparam.freq_per_poweplant,
        Aᵀ = Aᵀ,
        em = emᵀ
    )

end


"""
    update_c_pp(pp::PowerPlant, p_f::Float64)

Updates the marginal production cost of the power plant based on the type and 
    price of fossil fuels.
    Lamperti (2018), eq 11.
"""
function update_c_pp!(
    pp::PowerPlant, 
    p_f::Float64,
    τᶜ::Float64
    )

    pp.c = pp.type == "Dirty" ? p_f / pp.Aᵀ + pp.em * τᶜ : 0.0
end

function update_age_pp!(
    pp::PowerPlant
    )

    pp.age += 1
end