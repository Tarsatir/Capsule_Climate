"""
Climate box, containing:
    - time series on emissions
    - time series on carbon concentrations and net primary production
    - time series on temperature
    - parameters used to compute mechanics

Parameters and parameter values all follow from table 9 in Lamperti et al. (2018).
"""
@with_kw mutable struct Climate

    T::Int64

    # Emissions
    carbon_emissions::Vector{Float64} = zeros(Float64, T)      # carbon emissions
    carbon_emissions_kp::Vector{Float64} = zeros(Float64, T)   # carbon emissions of cp
    carbon_emissions_cp::Vector{Float64} = zeros(Float64, T)   # carbon emissions of kp
    emissions_index::Vector{Float64} = fill(100.0, T)          # index of carbon emissions
    energy_percentage::Vector{Float64} = zeros(Float64, T)     # percentage of emissions coming from energy
end


"""
Collect carbon emissions from producers.
"""
function collect_emissions_cl!(
    climate::Climate,
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    ep,
    t::Int64,
    t_warmup::Int64,
    model::ABM
    )

    climate.carbon_emissions_cp[t] = sum(cp_id -> model[cp_id].emissions, all_cp)
    climate.carbon_emissions_kp[t] = sum(kp_id -> model[kp_id].emissions, all_kp)
    climate.carbon_emissions[t] = (climate.carbon_emissions_kp[t] + 
                                   climate.carbon_emissions_cp[t] + 
                                   ep.emissions[t])

    # Start saving emission indeces from warmup time
    if t > t_warmup
        climate.emissions_index[t] = 100 * climate.carbon_emissions[t] / climate.carbon_emissions[t_warmup]
        climate.energy_percentage[t] = 1 - (climate.carbon_emissions_cp[t] + climate.carbon_emissions_kp[t]) / climate.carbon_emissions[t]
    end
end