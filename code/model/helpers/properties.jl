"""
Model properties
"""
mutable struct Properties

    t::Int64
    t_warmup::Int64
    T::Int64

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
    climate::Climate


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
    mdata_tosave::Union{Nothing, Vector{Symbol}}
    epdata_tosave::Union{Nothing, Vector{Symbol}}
    climatedata_tosave::Union{Nothing, Vector{Symbol}}
    governmentdata_tosave::Union{Nothing, Vector{Symbol}}

    macroeconomy::MacroEconomy
    labormarket::LaborMarket
    kp_brochures::Dict{Symbol, Dict{Symbol, Float64}}
    cmdata::CMData
    # ginidata::Matrix{Float64}
    perc_of_wealth::Vector{Float64}
    equal_div::Vector{Float64}
end