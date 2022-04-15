"""
Climate box, containing:
    - time series on emissions
    - time series on carbon concentrations and net primary production
    - time series on temperature
    - parameters used to compute mechanics

Parameters and parameter values all follow from table 9 in Lamperti et al. (2018).
"""
@with_kw mutable struct Climate

    T::Int

    # Emissions
    carbon_emissions::Vector{Float64} = zeros(Float64, T)      # carbon emissions
    carbon_emissions_kp::Vector{Float64} = zeros(Float64, T)   # carbon emissions of cp
    carbon_emissions_cp::Vector{Float64} = zeros(Float64, T)   # carbon emissions of kp

    # Carbon Cycle
    C_a::Vector{Float64} = fill(830.0, T)                      # atmospheric carbon levels
    C_m::Vector{Float64} = zeros(Float64, T)                   # mixed-layer oceanic carbon levels
    C_m_star::Vector{Float64} = zeros(Float64, T)              # mixed-layer reference carbon concentration
    C_d::Vector{Float64} = fill(10_010.0, T)                   # deep-layer oceanic carbon levels
    ΔC_md::Vector{Float64} = zeros(Float64, T)                 # net flux of carbon from mixed to deep layer
    NPP::Vector{Float64} = zeros(Float64, T)                   # Net primary production
    ξ::Vector{Float64} = zeros(Float64, T)                     # Revelle/buffer factor

    # Temperature
    T_a::Vector{Float64} = fill(14.8, T)                       # Temperature in atmosphere
    T_m::Vector{Float64} = zeros(Float64, T)                   # Temperature in mixed ocean layer
    T_d::Vector{Float64} = zeros(Float64, T)                   # Temperature in deep ocean layer

    # Pre-industrial reference levels
    T_0::Int = 14                    # Pre-industrial global mean surface temp
    C_md_0::Int = 10_237             # Pre-industrial carbon in ocean
    C_a_0::Int = 590                 # Pre-industrial carbon level
    NPP_0::Int = 85_177              # Pre-industrial NPP
    ξ_0::Float64 = 9.7               # Reference Revelle/buffer factor

    # Parameters
    β_C::Float64 = 1                 # Response of primary production to carbon concentration
    β_TC::Float64 = -0.01            # Sensitivity of carbon uptake to temp by land
    β_T::Float64 = 0.003             # Sensitivity of carbon uptake to temp

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
    cl::Climate,
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    ep::EnergyProducer,
    t::Int,
    model::ABM
    )

    cl.carbon_emissions_cp[t] = sum(cp_id -> model[cp_id].emissions, all_cp)
    cl.carbon_emissions_kp[t] = sum(kp_id -> model[kp_id].emissions, all_kp)
    cl.carbon_emissions[t] = (cl.carbon_emissions_kp[t] + 
                                          cl.carbon_emissions_cp[t] + 
                                          ep.emissions[t])
end


"""
Update carbon concentrations C_a, C_m and C_d.
"""
function compute_carbon_concentrations_cl!(
    cl::Climate,
    t::Int
    )

    # Compute NPP
    compute_NPP_cl!(cl, t)

    # Compute Revelle factor
    compute_Revelle_cl!(cl, t)

    # Compute reference carbon concentration
    compute_C_m_star_cl!(cl, t)

    # Compute carbon concentration in mixed layer
    compute_C_m_cl!(cl, t)

    # Compute flux of mixed to deep layer
    compute_ΔC_md_cl!(cl, t)
end


"""
Computes the net primary production (NPP).
    Lamperti et al (2018) eq 22.
"""
function compute_NPP_cl!(
    cl::Climate,
    t::Int
    )

    cl.NPP[t] = cl.NPP_0 * (1 + cl.β_C * log(cl.C_a[t] / cl.C_a_0)) * (1 - β_TC * cl.T_m[t-1])
end


"""
Computes the Revelle factor ξ.
    Lamperti et al (2018) eq 24.
"""
function compute_Revelle_cl!(
    cl::Climate,
    t::Int
    )

    cl.ξ[t] = cl.ξ_0 + cl.δ * log(cl.C_a[t-1] / cl.C_a_0)
end


"""
Computes reference carbon concentration in mixed layer C*ₘ.
    Lamperti et al (2018) eq 25.
"""
function compute_C_m_star_cl!(
    cl::Climate, 
    t::Int
    )

    cl.C_m_star[t] = cl.C_m[0] * (1 - cl.β_T * cl.T_m[t-1])
end


"""
Computes carbon concentration in mixed layer Cₘ.
    Lamperti et al (2018) eq 23.
"""
function compute_C_m_cl!(
    cl::Climate,
    t::Int
    )

    cl.C_m[t] = cl.C_m_star[t] * (cl.C_a[t] / cl.C_a_0) ^ (1 / cl.ξ[t])
end


"""
Computes net flux from mixed to deep layer ΔC_md
"""
function compute_ΔC_md_cl!(
    cl::Climate, 
    t::Int
    )

    cl.ΔC_md[t] = cl.k_eddy * ((cl.C_m[t-1] / cl.d_m) - (cl.C_d[t-1] / cl.d_d)) / mean(cl.d_d, cl.d_m)
end