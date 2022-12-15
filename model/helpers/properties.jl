"""
Model properties
"""
mutable struct Properties

    """*******************************************************
    INITIALIZING PARAMETERS
    *******************************************************"""
    i_param::InitParam


    """*******************************************************
    GLOBAL PARAMETERS
    *******************************************************"""
    g_param::GlobalParam


    """*******************************************************
    SINGLE AGENTS
    *******************************************************"""
    gov::Government
    ep::EnergyProducer
    idxf::IndexFund


    """*******************************************************
    SCHEDULES
    *******************************************************"""
    all_hh::Vector{Int64}
    all_cp::Vector{Int64}
    all_kp::Vector{Int64}
    all_p::Vector{Int64}

    
    """*******************************************************
    DATA STRUCTURES
    *******************************************************"""
    macroeconomy::MacroEconomy
    # firmdata::DataFrame
    # householddata::DataFrame
    kp_brochures::Dict{Symbol, Dict{Symbol, Float64}}
    cmdata::CMData

end