function compute_growthrates(
    ts::Vector{Float64}
    )

    return 100 .* (ts[1:end-1] .- ts[2:end]) ./ ts[2:end]
end


"""
Computes moments of runoutput and writes to dataframe or array
"""
function convertrunoutput(
    runoutput::RunOutput,
    sim_nr::Int64;
    return_as_df::Bool=false,
    t_warmup::Int64=100,
    )

    # Prepare data to be written to dataframe
    GDP_1st = mean(runoutput.GDP_growth[t_warmup:end])
    GDP_2nd = var(runoutput.GDP_growth[t_warmup:end])
    GDP_3rd = skewness(runoutput.GDP_growth[t_warmup:end])
    GDP_4th = kurtosis(runoutput.GDP_growth[t_warmup:end])

    acorr_GDP = cor(runoutput.GDP_growth[t_warmup+1:end], runoutput.GDP_growth[t_warmup:end-1])

    # Write unemployment data to dataframe
    U_1st = mean(runoutput.U[t_warmup:end])
    U_2nd = var(runoutput.U[t_warmup:end])

    dU1st = mean(runoutput.dU[t_warmup:end])
    corr_GDP_dU = cor(runoutput.GDP_growth[t_warmup:end], runoutput.dU[t_warmup:end])

    # Write Gini data to dataframe
    GINI_I_1st = mean(runoutput.GINI_I[t_warmup:end])
    GINI_I_2nd = var(runoutput.GINI_I[t_warmup:end])

    GINI_W_1st = mean(runoutput.GINI_W[t_warmup:end])
    GINI_W_2nd = var(runoutput.GINI_W[t_warmup:end])

    # Add productivity growth
    LP_g = compute_growthrates(runoutput.avg_π_LP)
    LP_g_1st = mean(LP_g[t_warmup:end])
    LP_g_2nd = var(LP_g[t_warmup:end])

    EE_g = compute_growthrates(runoutput.avg_π_EE)
    EE_g_1st = mean(EE_g[t_warmup:end])
    EE_g_2nd = var(EE_g[t_warmup:end])

    EF_g = compute_growthrates(runoutput.avg_π_EF)
    EF_g_1st = mean(EF_g[t_warmup:end])
    EF_g_2nd = var(EF_g[t_warmup:end])

    # Write poverty data to dataframe
    FGT_1st = mean(runoutput.FGT[t_warmup:end])
    FGT_2nd = var(runoutput.FGT[t_warmup:end])

    # Write emissions indexes
    em2030 = runoutput.emissions_index[220]
    em2040 = runoutput.emissions_index[340]
    em2050 = runoutput.emissions_index[460]

    if return_as_df
        return DataFrame(
                    :sim_nr => sim_nr,
                    :GDP_1st => GDP_1st,
                    :GDP_2nd => GDP_2nd,
                    :GDP_3rd => GDP_3rd,
                    :GDP_4th => GDP_4th,
                    :acorr_GDP => acorr_GDP,
                    :U_1st => U_1st,
                    :U_2nd => U_2nd,
                    :dU1st => dU1st,
                    :corr_GDP_dU => corr_GDP_dU,
                    :GINI_I_1st => GINI_I_1st,
                    :GINI_I_2nd => GINI_I_2nd,
                    :GINI_W_1st => GINI_W_1st,
                    :GINI_W_2nd => GINI_W_2nd,
                    :LP_g_1st => LP_g_1st,
                    :LP_g_2nd => LP_g_2nd,
                    :EE_g_1st => EE_g_1st,
                    :EE_g_2nd => EE_g_2nd,
                    :EF_g_1st => EF_g_1st,
                    :EF_g_2nd => EF_g_2nd,
                    :FGT_1st => FGT_1st,
                    :FGT_2nd => FGT_2nd,
                    :em2030 => em2030,
                    :em2040 => em2040,
                    :em2050 => em2050
                )
    else
        return [sim_nr,
                GDP_1st, GDP_2nd, GDP_3rd, GDP_4th, acorr_GDP,
                U_1st, U_2nd, dU1st, corr_GDP_dU, 
                GINI_I_1st, GINI_I_2nd, GINI_W_1st, GINI_W_2nd,
                LP_g_1st, LP_g_2nd, EE_g_1st, EE_g_2nd,
                EF_g_1st, EF_g_2nd, FGT_1st, FGT_2nd,
                em2030, em2040, em2050
                ]
    end
end