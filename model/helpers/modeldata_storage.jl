"""
SIMDATA - DATA STRUCT FOR DATA USED DURING SIMULATION
"""

@with_kw mutable struct SimData
    model::ABM 
    globalparam::GlobalParam 
    initparam::InitParam 
    macroeconomy::MacroEconomy 
    government::Government 
    ep::EnergyProducer 
    labormarket::LaborMarket 
    indexfund::IndexFund 
    climate::Climate 
    cmdata::CMData
    firmdata::Union{Nothing, DataFrame} = nothing
end


function gensimdata(
    model::ABM, 
    globalparam::GlobalParam, 
    initparam::InitParam,
    macroeconomy::MacroEconomy, 
    government::Government, 
    ep::EnergyProducer, 
    labormarket::LaborMarket, 
    indexfund::IndexFund, 
    climate::Climate, 
    cmdata::CMData;
    firmdata=nothing
    )

    return SimData(
        model,
        globalparam,
        initparam,
        macroeconomy,
        government,
        ep,
        labormarket,
        indexfund,
        climate,
        cmdata,
        firmdata
    )
end


"""
RUNOUTPUT - DATA STRUCT FOR FINAL AGGREGATE OUTPUT DATA
"""

"""
Struct that holds output of run when returned at end of simulation
"""
struct RunOutput
    GDP::Vector{Float64}
    GDP_growth::Vector{Float64}
    total_Q_growth::Vector{Float64}
    total_Q_cp::Vector{Float64}
    total_Q_kp::Vector{Float64}
    LIS::Vector{Float64}
    U::Vector{Float64}
    dU::Vector{Float64}
    C::Vector{Float64}
    I::Vector{Float64}
    wages::Vector{Float64}
    prices::Vector{Float64}
    markups::Vector{Float64}
    TotDebt::Vector{Float64}
    RD::Vector{Float64}
    EnDem::Vector{Float64}
    inventories::Vector{Float64}
    GINI_I::Vector{Float64}
    GINI_W::Vector{Float64}
    I_20::Vector{Float64}
    I_80::Vector{Float64}
    W_20::Vector{Float64}
    W_80::Vector{Float64}
    bankrupty_cp::Vector{Float64}
    avg_π_LP::Vector{Float64}
    avg_π_EE::Vector{Float64}
    avg_π_EF::Vector{Float64}
    emissions_total::Vector{Float64}
    emissions_index::Vector{Float64}
end

function genrunoutput(macroeconomy, ep, climate)
    return RunOutput(
        macroeconomy.GDP,
        macroeconomy.GDP_growth,
        macroeconomy.total_Q_growth,
        macroeconomy.total_Q_cp,
        macroeconomy.total_Q_kp,
        macroeconomy.LIS,
        macroeconomy.U,
        macroeconomy.dU,
        macroeconomy.total_C,
        macroeconomy.total_I,
        macroeconomy.w̄_avg,
        macroeconomy.p̄,
        macroeconomy.μ_cp,
        macroeconomy.debt_tot,
        macroeconomy.RD_total,
        ep.Dₑ,
        macroeconomy.N_goods,
        macroeconomy.GINI_I,
        macroeconomy.GINI_W,
        macroeconomy.I_20,
        macroeconomy.I_80,
        macroeconomy.W_20,
        macroeconomy.W_80,
        macroeconomy.bankrupt_cp,
        macroeconomy.avg_π_LP,
        macroeconomy.avg_π_EE,
        macroeconomy.avg_π_EF,
        climate.carbon_emissions,
        climate.emissions_index
    )
end
