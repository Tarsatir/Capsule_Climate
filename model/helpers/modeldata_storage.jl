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
    I_min::Vector{Float64}
    W_20::Vector{Float64}
    W_80::Vector{Float64}
    W_min::Vector{Float64}
    bankrupty_cp::Vector{Float64}
    avg_π_LP::Vector{Float64}
    avg_π_EE::Vector{Float64}
    avg_π_EF::Vector{Float64}
    Em::Vector{Float64}
    EmIndex::Vector{Float64}
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
        macroeconomy.I_min,
        macroeconomy.W_20,
        macroeconomy.W_80,
        macroeconomy.W_min,
        macroeconomy.bankrupt_cp,
        macroeconomy.avg_π_LP,
        macroeconomy.avg_π_EE,
        macroeconomy.avg_π_EF,
        climate.carbon_emissions,
        climate.emissions_index
    )
end


function genfirmdata(
    all_cp::Vector{Int64}, 
    all_kp::Vector{Int64}
    )::DataFrame

    start_cp = minimum(all_cp)
    end_cp = maximum(all_cp)
    start_kp = minimum(all_kp)
    end_kp = maximum(all_kp)

    # Create column names
    ages_cp = [Symbol("agecp_$i") for i in start_cp:end_cp]
    ages_kp = [Symbol("agekp_$i") for i in start_kp:end_kp]

    firmsize_cp = [Symbol("Scp_$i") for i in start_cp:end_cp]
    firmsize_kp = [Symbol("Skp_$i") for i in start_kp:end_kp]

    LP_cp = [Symbol("LPcp_$i") for i in start_cp:end_cp]
    LP_kp = [Symbol("LPkp_$i") for i in start_kp:end_kp]

    EE_cp = [Symbol("EEcp_$i") for i in start_cp:end_cp]
    EE_kp = [Symbol("EEkp_$i") for i in start_kp:end_kp]

    EF_cp = [Symbol("EFcp_$i") for i in start_cp:end_cp]
    EF_kp = [Symbol("EFkp_$i") for i in start_kp:end_kp]

    I_cp = [Symbol("Icp_$i") for i in start_cp:end_cp]
    # I_kp = [Symbol("Ikp_$i") for i in start_kp:end_kp]

    K_cp = [Symbol("Kcp_$i") for i in start_cp:end_cp]

    columnnames = vcat(
        [:t],
        ages_cp,
        ages_kp,
        firmsize_cp,
        firmsize_kp,
        LP_cp,
        LP_kp,
        EE_cp,
        EE_kp,
        EF_cp,
        EF_kp,
        I_cp,
        # I_kp
        K_cp
    )

    return DataFrame([[] for _ in columnnames], columnnames)
end


function appendfirmdata!(
    firmdata::DataFrame,
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    t::Int64,
    model::ABM
    )::DataFrame

    start_cp = minimum(all_cp)
    end_cp = maximum(all_cp)
    start_kp = minimum(all_kp)
    end_kp = maximum(all_kp)

    ages_cp = [model[cp_id].age for cp_id in start_cp:end_cp]
    ages_kp = [model[kp_id].age for kp_id in start_kp:end_kp]

    firmsize_cp = [model[cp_id].curracc.S for cp_id in start_cp:end_cp]
    firmsize_kp = [model[kp_id].curracc.S for kp_id in start_kp:end_kp]

    LP_cp = [model[cp_id].π_LP for cp_id in start_cp:end_cp]
    LP_kp = [model[kp_id].B_LP for kp_id in start_kp:end_kp]

    EE_cp = [model[cp_id].π_EE for cp_id in start_cp:end_cp]
    EE_kp = [model[kp_id].B_EE for kp_id in start_kp:end_kp]

    EF_cp = [model[cp_id].π_EF for cp_id in start_cp:end_cp]
    EF_kp = [model[kp_id].B_EF for kp_id in start_kp:end_kp]

    I_cp = [model[cp_id].curracc.TCI for cp_id in start_cp:end_cp]
    # I_kp = [model[cp_id].curracc.TCI for cp_id in start_kp:end_kp]

    K_cp = [model[cp_id].balance.K for cp_id in start_cp:end_cp]

    columnvalues = vcat(
        [t],
        ages_cp,
        ages_kp,
        firmsize_cp,
        firmsize_kp,
        LP_cp,
        LP_kp,
        EE_cp,
        EE_kp,
        EF_cp,
        EF_kp,
        I_cp,
        # I_kp
        K_cp
    )

    push!(firmdata, columnvalues)
    return firmdata
end