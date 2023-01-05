"""
File used to write simulation results to data files
"""

using DataFrames
using CSV

# """
# Saves macro variables of interest to csv

# Receives:
#     macroeconomy: mut struct with macro variables of interest
# """
# function save_macro_data(macroeconomy)

#     df = DataFrame(
#         GDP = macroeconomy.GDP,
#         GDP_I = macroeconomy.GDP_I,
#         GDP_cp = macroeconomy.GDP_Π_cp,
#         GDP_kp = macroeconomy.GDP_Π_kp,
#         GDP_growth = macroeconomy.GDP_growth,

#         total_C = macroeconomy.total_C,
#         total_C_actual = macroeconomy.total_C_actual,
#         total_I = macroeconomy.total_I,
#         total_w = macroeconomy.total_w,

#         LIS = macroeconomy.LIS,

#         returns_investments = macroeconomy.returns_investments,

#         unsat_demand = macroeconomy.unsat_demand,
#         unspend_C = macroeconomy.unspend_C,
#         unsat_invest = macroeconomy.unsat_invest,
#         unsat_L_demand = macroeconomy.unsat_L_demand,
#         avg_N_goods = macroeconomy.avg_N_goods,

#         p_avg_cp=macroeconomy.p̄,
#         CPI=macroeconomy.CPI,
#         CPI_kp = macroeconomy.CPI_kp,

#         mu_cp = macroeconomy.μ_cp,
#         mu_kp = macroeconomy.μ_kp,

#         M = macroeconomy.M,
#         M_hh = macroeconomy.M_hh,
#         M_cp = macroeconomy.M_cp,
#         M_kp = macroeconomy.M_kp,
#         M_ep = macroeconomy.M_ep,
#         M_gov = macroeconomy.M_gov,
#         M_if = macroeconomy.M_if,

#         debt_tot = macroeconomy.debt_tot,
#         debt_cp = macroeconomy.debt_cp,
#         debt_cp_allowed = macroeconomy.debt_cp_allowed,
#         debt_kp = macroeconomy.debt_kp,
#         debt_kp_allowed = macroeconomy.debt_kp_allowed,
#         debt_unpaid_kp = macroeconomy.debt_unpaid_kp,
#         debt_unpaid_cp = macroeconomy.debt_unpaid_cp,

#         UR = macroeconomy.U,
#         switch_rate = macroeconomy.switch_rate,
#         Exp_UB=macroeconomy.Exp_UB,

#         s_emp = macroeconomy.s̄_emp,
#         s_unemp = macroeconomy.s̄_unemp,

#         w_avg = macroeconomy.w̄_avg,
#         wr_avg = macroeconomy.wʳ_avg,
#         ws_avg = macroeconomy.wˢ_avg,
#         wo_max_avg = macroeconomy.wᴼ_max_mean,

#         I_avg = macroeconomy.Ī_avg,
#         I_labor_avg = macroeconomy.I_labor_avg,
#         I_capital_avg = macroeconomy.I_capital_avg,
#         I_UB_avg = macroeconomy.I_UB_avg,
#         I_socben_avg = macroeconomy.I_socben_avg,

#         dL_avg = macroeconomy.ΔL̄_avg,
#         dL_std = macroeconomy.ΔL̄_std,
#         dL_cp_avg = macroeconomy.ΔL̄_cp_avg,
#         dL_kp_avg = macroeconomy.ΔL̄_kp_avg,

#         EI_avg = macroeconomy.EI_avg,
#         n_mach_EI = macroeconomy.n_mach_EI_avg,
#         RS_avg = macroeconomy.RS_avg,
#         n_mach_RS = macroeconomy.n_mach_RS_avg,

#         avg_pi_LP = macroeconomy.avg_π_LP,
#         avg_pi_EE = macroeconomy.avg_π_EE,
#         avg_pi_EF = macroeconomy.avg_π_EF,

#         avg_A_LP = macroeconomy.avg_A_LP,
#         avg_A_EE = macroeconomy.avg_A_EE,
#         avg_A_EF = macroeconomy.avg_A_EF,

#         avg_B_LP = macroeconomy.avg_B_LP,
#         avg_B_EE = macroeconomy.avg_B_EE,
#         avg_B_EF = macroeconomy.avg_B_EF,

#         total_Q_cp = macroeconomy.total_Q_cp,
#         total_Q_kp = macroeconomy.total_Q_kp,

#         avg_Q_cp = macroeconomy.avg_Q_cp,
#         avg_Qs_cp = macroeconomy.avg_Qˢ_cp,
#         avg_Qe_cp = macroeconomy.avg_Qᵉ_cp,
#         avg_Q_kp = macroeconomy.avg_Q_kp,
#         avg_D_cp = macroeconomy.avg_D_cp,
#         avg_Du_cp = macroeconomy.avg_Dᵁ_cp,
#         avg_De_cp = macroeconomy.avg_Dᵉ_cp,

#         bankrupt_cp = macroeconomy.bankrupt_cp,
#         bankrupt_kp = macroeconomy.bankrupt_kp,

#         cu = macroeconomy.cu,
#         avg_n_machines_cp = macroeconomy.avg_n_machines_cp,

#         gini_I = macroeconomy.GINI_I,
#         gini_W = macroeconomy.GINI_W,

#         I_min = macroeconomy.I_min,
#         I_20 = macroeconomy.I_20,
#         I_80 = macroeconomy.I_80,
#         I_max = macroeconomy.I_max,

#         W_min = macroeconomy.W_min,
#         W_20 = macroeconomy.W_20,
#         W_80 = macroeconomy.W_80,
#         W_max = macroeconomy.W_max
#     )
#     CSV.write("results/result_data/first.csv", df)

#     CSV.write("results/result_data/alpha_W_quantiles.csv", DataFrame(macroeconomy.α_W_quantiles, :auto))
# end


function save_simdata(
    agent_df::DataFrame,
    model_df::DataFrame,
    seed::Int64,
)
    # NOTE: CONVERSION DF TO STRING IS TMP SOLUTION, SHOULD BE FIXED BACK WHEN PACKAGES 
        #       ARE CONSISTENT AGAIN!
    CSV.write(string("results/result_data/agent_data_", seed, ".csv"), string.(agent_df))
    CSV.write(string("results/result_data/model_data_", seed, ".csv"), model_df)
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

        p_e = energy_producer.p_ep,

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
    )
    CSV.write("results/result_data/climate_and_energy.csv", df)
end


function save_household_quartiles(
    householddata::Array
)

    CSV.write("results/result_data/household_quantiles.csv", householddata[2])
end