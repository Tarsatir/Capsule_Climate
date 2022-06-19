"""
File used to write simulation results to data files
"""

using DataFrames
using CSV

"""
Saves macro variables of interest to csv

Receives:
    macro_struct: mut struct with macro variables of interest
"""
function save_macro_data(macro_struct)

    df = DataFrame(
        GDP = macro_struct.GDP,
        GDP_I = macro_struct.GDP_I,
        GDP_cp = macro_struct.GDP_Π_cp,
        GDP_kp = macro_struct.GDP_Π_kp,
        GDP_growth = macro_struct.GDP_growth,

        C = macro_struct.C,

        unsat_demand = macro_struct.unsat_demand,
        unspend_C = macro_struct.unspend_C,
        unsat_invest = macro_struct.unsat_invest,
        unsat_L_demand = macro_struct.unsat_L_demand,
        avg_N_goods = macro_struct.avg_N_goods,

        CPI=macro_struct.CPI,
        CPI_kp = macro_struct.CPI_kp,

        # mu_bp = macro_struct.μ_bp,
        # mu_lp = macro_struct.μ_lp,
        mu_cp = macro_struct.μ_cp,
        mu_kp = macro_struct.μ_kp,

        M = macro_struct.M,
        M_hh = macro_struct.M_hh,
        M_cp = macro_struct.M_cp,
        M_kp = macro_struct.M_kp,
        M_ep = macro_struct.M_ep,
        M_gov = macro_struct.M_gov,
        M_if = macro_struct.M_if,

        debt_tot = macro_struct.debt_tot,
        debt_cp = macro_struct.debt_cp,
        debt_cp_allowed = macro_struct.debt_cp_allowed,
        debt_kp = macro_struct.debt_kp,
        debt_kp_allowed = macro_struct.debt_kp_allowed,
        debt_unpaid_kp = macro_struct.debt_unpaid_kp,
        debt_unpaid_cp = macro_struct.debt_unpaid_cp,

        UR = macro_struct.U,
        switch_rate = macro_struct.switch_rate,
        Exp_UB=macro_struct.Exp_UB,

        s_emp = macro_struct.s̄_emp,
        s_unemp = macro_struct.s̄_unemp,

        w_avg = macro_struct.w̄_avg,
        wr_avg = macro_struct.wʳ_avg,
        ws_avg = macro_struct.wˢ_avg,
        wo_max_avg = macro_struct.wᴼ_max_mean,

        I_avg = macro_struct.Ī_avg,
        I_labor_avg = macro_struct.I_labor_avg,
        I_capital_avg = macro_struct.I_capital_avg,
        I_UB_avg = macro_struct.I_UB_avg,
        I_socben_avg = macro_struct.I_socben_avg,

        dL_avg = macro_struct.ΔL̄_avg,
        dL_std = macro_struct.ΔL̄_std,
        dL_cp_avg = macro_struct.ΔL̄_cp_avg,
        dL_kp_avg = macro_struct.ΔL̄_kp_avg,

        EI_avg = macro_struct.EI_avg,
        n_mach_EI = macro_struct.n_mach_EI_avg,
        RS_avg = macro_struct.RS_avg,
        n_mach_RS = macro_struct.n_mach_RS_avg,

        avg_pi_LP =  macro_struct.avg_π_LP,
        avg_pi_EE = macro_struct.avg_π_EE,

        avg_A_LP = macro_struct.avg_A_LP,
        avg_A_EE = macro_struct.avg_A_EE,
        avg_A_EF = macro_struct.avg_A_EF,

        avg_B_LP = macro_struct.avg_B_LP,
        avg_B_EE = macro_struct.avg_B_EE,
        avg_B_EF = macro_struct.avg_B_EF,

        avg_Q_cp = macro_struct.avg_Q_cp,
        avg_Qs_cp = macro_struct.avg_Qˢ_cp,
        avg_Qe_cp = macro_struct.avg_Qᵉ_cp,
        avg_Q_kp = macro_struct.avg_Q_kp,
        avg_D_cp = macro_struct.avg_D_cp,
        avg_Du_cp = macro_struct.avg_Dᵁ_cp,
        avg_De_cp = macro_struct.avg_Dᵉ_cp,

        bankrupt_cp = macro_struct.bankrupt_cp,
        bankrupt_kp = macro_struct.bankrupt_kp,

        cu = macro_struct.cu,
        avg_n_machines_cp = macro_struct.avg_n_machines_cp,

        gini_I = macro_struct.GINI_I,
        gini_W = macro_struct.GINI_W,
        FGT = macro_struct.FGT
    )
    CSV.write("results/result_data/first.csv", df)
end


function save_final_dist(
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    all_kp::Vector{Int}, 
    model::ABM
    )

    # Save income data of households 
    df = DataFrame(
        all_I = map(hh_id -> model[hh_id].total_I, all_hh),
        all_w = map(hh_id -> model[hh_id].w[end], all_hh),
        all_W = map(hh_id -> model[hh_id].W, all_hh),
        skills = map(hh_id -> model[hh_id].skill, all_hh)
    )
    CSV.write("results/result_data/final_income_dists.csv", df)

    # Save sales, profits and market share of cp
    df = DataFrame(
        all_S_cp = map(cp_id -> model[cp_id].curracc.S, all_cp),
        all_profit_cp = map(cp_id -> model[cp_id].Π[end], all_cp),
        all_f_cp = map(cp_id -> model[cp_id].f[end], all_cp),
        all_L_cp = map(cp_id -> model[cp_id].L, all_cp),
        all_p_cp = map(cp_id -> model[cp_id].p[end], all_cp),
        all_w_cp = map(cp_id -> model[cp_id].w̄[end], all_cp)
    )
    CSV.write("results/result_data/final_profit_dists_cp.csv", df)

    # Save sales, profits and market share of kp
    df = DataFrame(
        all_S_kp = map(kp_id -> model[kp_id].curracc.S, all_kp),
        all_profit_kp = map(kp_id -> model[kp_id].Π[end], all_kp),
        all_f_kp = map(kp_id -> model[kp_id].f[end], all_kp),
        all_L_kp = map(kp_id -> model[kp_id].L, all_kp)
    )
    CSV.write("results/result_data/final_profit_dists_kp.csv", df)

end


function save_climate_data(
    energy_producer,
    climate,
    model::ABM
    )

    df = DataFrame(
        energy_demand = energy_producer.Dₑ,
        total_capacity = energy_producer.Q̄ₑ,
        green_capacity = energy_producer.green_capacity,
        dirty_capacity = energy_producer.dirty_capacity,

        p_e = energy_producer.pₑ,

        RD = energy_producer.RDₑ,
        IN_g = energy_producer.IN_g,
        IN_d = energy_producer.IN_d,

        IC_g = energy_producer.IC_g,
        A_d = energy_producer.Aᵀ_d,
        em_d = energy_producer.emᵀ_d,
        c_d = energy_producer.c_d,

        emissions_total = climate.carbon_emissions,
        emissions_kp = climate.carbon_emissions_kp,
        emissions_cp = climate.carbon_emissions_cp,
        emissions_ep = energy_producer.emissions,

        C_a = climate.C_a,
        C_m = climate.C_m,
        C_d = climate.C_d,

        NPP = climate.NPP,

        dT_m = climate.δT_m,
        dT_d = climate.δT_d
    )
    CSV.write("results/result_data/climate_and_energy.csv", df)
end