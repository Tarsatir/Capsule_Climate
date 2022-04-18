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
    ΔC_am::Vector{Float64} = zeros(Float64, T)                 # net flux of carbon from atmosphere to mixed layer
    ΔC_md::Vector{Float64} = zeros(Float64, T)                 # net flux of carbon from mixed to deep layer
    NPP::Vector{Float64} = zeros(Float64, T)                   # Net primary production
    ξ::Vector{Float64} = zeros(Float64, T)                     # Revelle/buffer factor

    # Temperature
    F_CO2::Vector{Float64} = zeros(Float64, T)                 # Radiative forcing over time
    δT_a::Vector{Float64} = fill(14.8, T)                       # Temperature anomaly in atmosphere wrt pre-industrial
    δT_m::Vector{Float64} = zeros(Float64, T)                   # Temperature anomaly in mixed ocean layer wrt pre-industrial
    δT_d::Vector{Float64} = zeros(Float64, T)                   # Temperature anomaly in deep ocean layer wrt pre-industrial

    # Pre-industrial reference levels
    T_m_0::Int = 14                  # Pre-industrial global mean surface temp
    T_d_0::Int = 4                   # Pre-industrial deep ocean temp                 
    C_md_0::Int = 10_237e9           # Pre-industrial carbon in ocean
    C_a_0::Int = 590e9               # Pre-industrial carbon level
    NPP_0::Float64 = 85.177e9 / 4    # Pre-industrial NPP
    ξ_0::Float64 = 9.7               # Reference Revelle/buffer factor

    # Initial levels
    δT_a_init::Float64 = 14.8 - T_m_0 # Initial global mean surface temp anaomly wrt pre-industrial
    δT_m_init::Float64 = 14.8 - T_m_0 # Initial global mean surface temp anaomly wrt pre-industrial
    δT_d_init::Float64 = 4.8 - T_d_0  # Initial deep ocean layer temp anaomly wrt pre-industrial
    C_a_init::Float64 = 830e9         # Initial carbon concentration in atmosphere
    C_m_init::Float64 = 830e9         # Initial carbon concentration in mixed layer
    C_d_init::Float64 = 10_010e9      # Initial carbon concentration in deep layer     

    # Parameters
    β_C::Float64 = 1                 # Response of primary production to carbon concentration
    β_TC::Float64 = -0.01            # Sensitivity of carbon uptake to temp by land
    β_T::Float64 = 0.003             # Sensitivity of carbon uptake to temp

    δ::Float64 = 3.92                # Index response of buffer factor to carbon concentration
    σ::Float64 = 1.0                 # transfer rate of water from upper to lower ocean layer

    k_eddy::Int = 1                  # Eddy diffusion coefficient
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
    ep,
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
Update carbon concentrations C_a, C_m and C_d, updates temperature.
    Based on section 2.3.3 in Lamperti et al (2018)
"""
function carbon_equilibrium_tempchange_cl!(
    cl::Climate,
    t::Int
    )

    g = 1e5

    # (1) Determine the atmospheric carbon amount without interactions with the
    # biospehere or oceans
    C_a_t = t > 1 ? cl.C_a[t-1] + cl.carbon_emissions[t] * g : cl.C_a_init + cl.carbon_emissions[t] * g

    # (2) Determine new capacity of oceans to update CO2

    # Compute Revelle factor
    compute_Revelle_cl!(cl, t)

    # Compute reference carbon concentration
    compute_C_m_star_cl!(cl, t)

    # (3) Carbon exchanges between atmosphere and biosphere and oceans

    # Compute NPP
    compute_NPP_cl!(cl, C_a_t, t)

    # Compute flux of mixed to deep layer
    compute_ΔC_md_cl!(cl, t)

    # Compute carbon concentration in mixed layer
    compute_C_m_cl!(cl, C_a_t, t)

    # (4) Compute the new carbon equilibria
    compute_carbon_equilibria_cl!(cl, C_a_t, t)

    # (5)-(6) Compute new radiative forcing and new temperature
    compute_F_CO2_cl!(cl, t)
    compute_new_temp_cl!(cl, t)
end


"""
Computes the net primary production (NPP).
    Lamperti et al (2018) eq 22.
"""
function compute_NPP_cl!(
    cl::Climate,
    C_a_t::Float64,
    t::Int
    )

    # println(C_a_t)
    # println(cl.C_a_init)
    # println(log(C_a_t / cl.C_a_init))

    if t > 1
        cl.NPP[t] = cl.NPP_0 * (1 + cl.β_C * log(C_a_t / cl.C_a_init)) * (1 - cl.β_TC * cl.δT_m[t-1])
    else
        cl.NPP[t] = cl.NPP_0 * (1 + cl.β_C * log(C_a_t / cl.C_a_init)) * (1 - cl.β_TC * cl.δT_m_init)
    end
end


"""
Computes the Revelle factor ξ.
    Lamperti et al (2018) eq 24.
"""
function compute_Revelle_cl!(
    cl::Climate,
    t::Int
    )

    cl.ξ[t] = t > 1 ? cl.ξ_0 + cl.δ * log(cl.C_a[t-1] / cl.C_a_init) : cl.ξ_0
end


"""
Computes reference carbon concentration in mixed layer C*ₘ.
    Lamperti et al (2018) eq 25.
"""
function compute_C_m_star_cl!(
    cl::Climate, 
    t::Int
    )

    cl.C_m_star[t] = t > 1 ? cl.C_m_init * (1 - cl.β_T * cl.δT_m[t-1]) : cl.C_m_init
end


"""
Computes carbon concentration in mixed layer Cₘ.
    Lamperti et al (2018) eq 23.
"""
function compute_C_m_cl!(
    cl::Climate,
    C_a_t::Float64,
    t::Int
    )

    cl.C_m[t] = cl.C_m_star[t] * (C_a_t / cl.C_a_init) ^ (1 / cl.ξ[t])
    cl.ΔC_am[t] = t > 1 ? cl.C_m[t] - cl.C_m[t-1] : cl.C_m[t] - cl.C_m_init
end


"""
Computes net flux from mixed to deep layer ΔC_md.
    Lamperti et al (2018) eq 26.
"""
function compute_ΔC_md_cl!(
    cl::Climate, 
    t::Int
    )

    if t > 1
        cl.ΔC_md[t] = cl.k_eddy * (cl.C_m[t-1] / cl.d_mixed - cl.C_d[t-1] / cl.d_deep) / mean((cl.d_deep, cl.d_mixed))
    else
        cl.ΔC_md[t] = cl.k_eddy * (cl.C_m_init / cl.d_mixed - cl.C_d_init / cl.d_deep) / mean((cl.d_deep, cl.d_mixed))
    end
end


"""
Computes the new carbon dioxide equilibria in the atmosphere, mixed layer and deep layer.
"""
function compute_carbon_equilibria_cl!(
    cl::Climate, 
    C_a_t::Float64,
    t::Int
    )

    # Compute how much was taken up by biospehere and oceans
    carbon_uptaken = cl.NPP[t] + cl.ΔC_am[t]
    # carbon_uptaken = 0
    cl.C_a[t] = t > 1 ? C_a_t - carbon_uptaken : C_a_t - carbon_uptaken

    # Exchange carbon between mixed and deep layers
    cl.C_m[t] -= cl.ΔC_md[t]
    cl.C_d[t] = t > 1 ? cl.C_d[t-1] + cl.ΔC_md[t] : cl.C_d_init + cl.ΔC_md[t]
end


"""
Computes the new value for ratiate forcing F_CO2.
    Lamperti et al (2018) eq 29.
"""
function compute_F_CO2_cl!(
    cl::Climate, 
    t::Int
    )

    cl.F_CO2[t] = cl.γ * log(cl.C_a[t] / cl.C_a_init)
end


"""
Computes new equilibrium temperatures for atmosphere, mixed and deep ocean layers
    Lamperti et al (2018) eq 27 and 28.
"""
function compute_new_temp_cl!(
    cl::Climate, 
    t::Int
    )

    if t > 1
        cl.δT_m[t] = cl.δT_m[t-1] + cl.c1 * (cl.F_CO2[t] - cl.λ * cl.δT_m[t-1] - 
                    cl.c3 * (cl.δT_m[t-1] - cl.δT_d[t-1]))
        cl.δT_d[t] = cl.δT_d[t-1] + cl.c4 * (cl.σ * (cl.δT_m[t-1] - cl.δT_d[t-1]))
    else
        cl.δT_m[t] = cl.δT_m_init + cl.c1 * (cl.F_CO2[t] - cl.λ * cl.δT_m_init - 
                    cl.c3 * (cl.δT_m_init - cl.δT_d_init))
        cl.δT_d[t] = cl.δT_d_init + cl.c4 * (cl.σ * (cl.δT_m_init - cl.δT_d_init))
    end
end