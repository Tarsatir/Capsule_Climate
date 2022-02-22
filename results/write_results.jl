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
        M=macro_struct.M,
        M_hh=macro_struct.M_hh,
        M_cp=macro_struct.M_cp,
        M_kp=macro_struct.M_kp,
        M_gov=macro_struct.M_gov,

        Deb_tot=macro_struct.Deb_tot,
        Deb_cp=macro_struct.Deb_cp,
        Deb_kp=macro_struct.Deb_kp,

        UR=macro_struct.U,
        Exp_UB=macro_struct.Exp_UB,
        s_avg=macro_struct.s̄_avg,
        s_std=macro_struct.s̄_std,
        w_avg=macro_struct.w̄_avg,
        w_std=macro_struct.w̄_std,
        wr_avg=macro_struct.wʳ_avg,
        wr_std=macro_struct.wʳ_std,
        ws_avg=macro_struct.wˢ_avg,
        ws_std=macro_struct.wˢ_std,
        I_avg=macro_struct.Ī_avg,
        I_std=macro_struct.Ī_std,

        dL_avg=macro_struct.ΔL̄_avg,
        dL_std=macro_struct.ΔL̄_std,
        dL_cp_avg=macro_struct.ΔL̄_cp_avg,
        dL_kp_avg=macro_struct.ΔL̄_kp_avg,

        EI_avg=macro_struct.EI_avg,
        RS_avg=macro_struct.RS_avg,

        avg_pi=macro_struct.avg_π,
        avg_A=macro_struct.avg_A,
        avg_B=macro_struct.avg_B,

        avg_Q_bp=macro_struct.avg_Q_bp,
        avg_Q_lp=macro_struct.avg_Q_lp,
        avg_Q_kp=macro_struct.avg_Q_kp
        )
    CSV.write("results/result_data/first.csv", df)
end

function save_final_dist(all_hh, model)

    df = DataFrame(all_I = map(hh_id -> model[hh_id].I[end], all_hh),
                   all_w = map(hh_id -> model[hh_id].w[end], all_hh),
                   all_W = map(hh_id -> model[hh_id].W[end], all_hh))
    CSV.write("results/result_data/final_dists.csv", df)

end