"""
Climate box, containing:
    - time series on emissions
    - time series on carbon concentrations and net primary production
    - time series on temperature
    - parameters used to compute mechanics

Parameters and parameter values all follow from table 9 in Lamperti et al. (2018).
"""
@with_kw mutable struct Climate{V<:Vector{Float64}}

    T::Int64

    # Emissions
    carbon_emissions::V = zeros(Float64, T)      # carbon emissions
    carbon_emissions_cp::V = zeros(Float64, T)   # carbon emissions of cp
    carbon_emissions_kp::V = zeros(Float64, T)   # carbon emissions of kp
    carbon_emissions_ep::V = zeros(Float64, T)   # carbon emissions of ep

    em_index::V = fill(100., T)
    em_index_cp::V = fill(100., T)
    em_index_kp::V = fill(100., T)
    em_index_ep::V = fill(100., T)
    energy_percentage::V = zeros(Float64, T)     # percentage of emissions coming from energy
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
    climate.carbon_emissions_ep[t] = ep.emissions[t]
    climate.carbon_emissions[t] = (climate.carbon_emissions_cp[t] + 
                                   climate.carbon_emissions_kp[t] + 
                                   climate.carbon_emissions_ep[t])

    # Start saving emission indeces from warmup time
    # if t > t_warmup
        # climate.emissions_index[t] = 100 * climate.carbon_emissions[t] / climate.carbon_emissions[t_warmup]
        
    # end
    climate.energy_percentage[t] = climate.carbon_emissions_ep[t] / climate.carbon_emissions[t]
end

compute_index(ts, t_index) = (ts ./ ts[t_index]) .* 100

function compute_emission_indices!(
    climate::Climate,
    t_index::Int64
)
    climate.em_index .= compute_index(climate.carbon_emissions, t_index)
    climate.em_index_cp .= compute_index(climate.carbon_emissions_cp, t_index)
    climate.em_index_kp .= compute_index(climate.carbon_emissions_kp, t_index)
    climate.em_index_ep .= compute_index(climate.carbon_emissions_ep, t_index)
end
