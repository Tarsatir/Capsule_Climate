mutable struct FirmTimeSeries
    cp_data::DataFrame
    kp_data::DataFrame
end


function update_firm_time_series!(
    t::Int, 
    model::ABM,
    firm_time_series::FirmTimeSeries
)
    # Get current firm data
    current_cp_data = get_cp_mdata(model)
    current_kp_data = get_kp_mdata(model)

    # Initialize if empty
    if isempty(firm_time_series.cp_data)
        firm_time_series.cp_data = current_cp_data
        firm_time_series.kp_data = current_kp_data
        return
    end

    # Check if both DataFrames have the same column names for cp_data
    if names(firm_time_series.cp_data) != names(current_cp_data)
        error("Column names for cp_data do not match.")
    end

    # Check if both DataFrames have the same column names for kp_data
    if names(firm_time_series.kp_data) != names(current_kp_data)
        error("Column names for kp_data do not match.")
    end

    # Append to the time series
    append!(firm_time_series.cp_data, current_cp_data)
    append!(firm_time_series.kp_data, current_kp_data)
end


