@with_kw mutable struct Climate

    T::Int
    
    # Carbon Cycle
    carbon_emissions::Vector{Float64} = zeros(Float64, T)      # carbon emissions
    carbon_emissions_kp::Vector{Float64} = zeros(Float64, T)   # carbon emissions of cp
    carbon_emissions_cp::Vector{Float64} = zeros(Float64, T)   # carbon emissions of kp

    C_a::Vector{Float64} = fill(830.0, T)                      # atmospheric carbon levels
    C_m::Vector{Float64} = zeros(Float64, T)                   # mixed-layer oceanic carbon levels
    C_d::Vector{Float64} = fill(10_010.0, T)                   # deep-layer oceanic carbon levels

    # Temperature
    T_a::Vector{Float64} = fill(14.8, T)                       # Temperature in atmosphere


    # Pre-industrial reference levels
    T_0::Int = 14                    # Pre-industrial global mean surface temp
    C_md_0::Int = 10_237             # Pre-industrial carbon in ocean
    C_a_0::Int = 590                 # Pre-industrial carbon level
    NPP_0::Int = 85_177              # Pre-industrial NPP

    # Parameters, all following from table 9 in Lamperti et al. (2018)
    β_C::Float64 = 1                 # Response of primary production to carbon concentration
    β_TC::Float64 = -0.01            # Sensitivity of carbon uptake to temp by land
    β_T::Float64 = 0.003             # Sensitivity of carbon uptake to temp

    ξ::Float64 = 9.7                 # Reference buffer factor
    δ::Float64 = 3.92                # Index response of buffer factor to carbon concentration

    d_eddy::Int = 1                  # Eddy diffusion coefficient
    d_mixed::Int = 100               # Mixed ocean depth
    d_deep::Int = 3500               # Deep ocean depth

    c1::Float64 = 0.098              # Diffusion for atmospheric temp eq
    c3::Float64 = 0.088              # Diffusion for deep oceans temp eq
    c4::Float64 = 0.025              # Sensitivity of atmospheric temp to deep ocean temp
    
    λ::Float64 = 2.9                 # Equilibrium climate sensitivity
    γ::Float64 = 5.35                # Radiative forcing coefficient


end


"""
Collect carbon emissions from producers.
"""
function collect_emissions_cl!(
    climate_struct::Climate,
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    ep::EnergyProducer,
    t::Int,
    model::ABM
    )

    climate_struct.carbon_emissions_cp[t] = sum(cp_id -> model[cp_id].emissions, all_cp)
    climate_struct.carbon_emissions_kp[t] = sum(kp_id -> model[kp_id].emissions, all_kp)
    climate_struct.carbon_emissions[t] = (climate_struct.carbon_emissions_kp[t] + 
                                          climate_struct.carbon_emissions_cp[t] + 
                                          ep.emissions[t])
end


"""
Update carbon concentrations C_a, C_m and C_d.
"""
function compute_carbon_concentrations_cl!(
    climate_struct::Climate,
    t::Int
    )


end