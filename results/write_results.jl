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
        GDP=macro_struct.GDP,
        GDP_I=macro_struct.GDP_I,
        GDP_cp=macro_struct.GDP_Π_cp,
        GDP_kp=macro_struct.GDP_Π_kp,
        GDP_growth=macro_struct.GDP_growth,

        C = macro_struct.C,

        unsat_demand = macro_struct.unsat_demand,
        avg_N_goods = macro_struct.avg_N_goods,

        CPI=macro_struct.CPI,
        CPI_kp = macro_struct.CPI_kp,

        mu_bp = macro_struct.μ_bp,
        mu_lp = macro_struct.μ_lp,
        mu_kp = macro_struct.μ_kp,

        M=macro_struct.M,
        M_hh=macro_struct.M_hh,
        M_cp=macro_struct.M_cp,
        M_kp=macro_struct.M_kp,
        M_gov=macro_struct.M_gov,
        M_if=macro_struct.M_if,

        debt_tot=macro_struct.debt_tot,
        debt_cp=macro_struct.debt_cp,
        debt_cp_allowed=macro_struct.debt_cp_allowed,
        debt_kp=macro_struct.debt_kp,
        debt_kp_allowed=macro_struct.debt_kp_allowed,
        debt_unpaid_kp = macro_struct.debt_unpaid_kp,
        debt_unpaid_cp = macro_struct.debt_unpaid_cp,

        UR=macro_struct.U,
        Exp_UB=macro_struct.Exp_UB,
        s_avg=macro_struct.s̄_avg,
        s_std=macro_struct.s̄_std,
        w_avg=macro_struct.w̄_avg,
        # w_std=macro_struct.w̄_std,
        wr_avg=macro_struct.wʳ_avg,
        # wr_std=macro_struct.wʳ_std,
        ws_avg=macro_struct.wˢ_avg,
        # ws_std=macro_struct.wˢ_std,
        wo_max_avg = macro_struct.wᴼ_max_mean,

        I_avg=macro_struct.Ī_avg,
        I_std=macro_struct.Ī_std,

        dL_avg=macro_struct.ΔL̄_avg,
        dL_std=macro_struct.ΔL̄_std,
        dL_cp_avg=macro_struct.ΔL̄_cp_avg,
        dL_kp_avg=macro_struct.ΔL̄_kp_avg,

        EI_avg=macro_struct.EI_avg,
        n_mach_EI=macro_struct.n_mach_EI_avg,
        RS_avg=macro_struct.RS_avg,
        n_mach_RS=macro_struct.n_mach_RS_avg,

        avg_pi=macro_struct.avg_π,
        avg_A=macro_struct.avg_A,
        avg_B=macro_struct.avg_B,

        avg_Q_bp=macro_struct.avg_Q_bp,
        avg_Q_lp=macro_struct.avg_Q_lp,
        avg_Q_kp=macro_struct.avg_Q_kp,

        bankrupt_bp = macro_struct.bankrupt_bp,
        bankrupt_lp = macro_struct.bankrupt_lp,
        bankrupt_kp = macro_struct.bankrupt_kp,

        cu = macro_struct.cu,
        avg_n_machines_bp = macro_struct.avg_n_machines_bp,
        avg_n_machines_lp = macro_struct.avg_n_machines_lp,

        gini_I = macro_struct.GINI_I,
        gini_W = macro_struct.GINI_W
    )
    CSV.write("results/result_data/first.csv", df)
end

function save_final_dist(
    all_hh::Vector{Int},
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    all_kp::Vector{Int}, 
    model::ABM
    )

    df = DataFrame(
            all_I = map(hh_id -> model[hh_id].I[end], all_hh),
            all_w = map(hh_id -> model[hh_id].w[end], all_hh),
            all_W = map(hh_id -> model[hh_id].W[end], all_hh),
    )
    CSV.write("results/result_data/final_income_dists.csv", df)

    df = DataFrame(
        all_S_bp = map(bp_id -> model[bp_id].curracc.S, all_bp),
        all_profit_bp = map(bp_id -> model[bp_id].Π[end], all_bp),
        all_S_lp = map(lp_id -> model[lp_id].curracc.S, all_lp),
        all_profit_lp = map(lp_id -> model[lp_id].Π[end], all_lp),
        all_f_bp = map(bp_id -> model[bp_id].f[end], all_bp),
        all_f_lp = map(lp_id -> model[lp_id].f[end], all_lp),
    )
    CSV.write("results/result_data/final_profit_dists_cp.csv", df)

    df = DataFrame(
        all_S_kp = map(kp_id -> model[kp_id].curracc.S, all_kp),
        all_profit_kp = map(kp_id -> model[kp_id].Π[end], all_kp),
        all_f_kp = map(kp_id -> model[kp_id].f[end], all_kp)
    )
    CSV.write("results/result_data/final_profit_dists_kp.csv", df)

end