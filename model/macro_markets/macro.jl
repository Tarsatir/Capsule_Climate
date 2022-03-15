
mutable struct MacroEconomy
    GDP :: Vector{Float64}       # GDP over time
    GDP_I :: Vector{Float64}     # income share of GDP over time
    GDP_Π_cp :: Vector{Float64}  # profit share of GDP of cp over time
    GDP_Π_kp :: Vector{Float64}  # profit share of GDP of kp over time
    p̄::Vector{Float64}           # average price over time
    CPI :: Vector{Float64}       # inflation over time
    C :: Vector{Float64}         # aggregate consumption over time

    # Division of money over sectors
    M :: Vector{Float64}         # total amount of money (should be stable)
    M_hh :: Vector{Float64}      # total amount of money at hh
    M_cp :: Vector{Float64}      # total amount of money at cp
    M_kp :: Vector{Float64}      # total amount of money at kp
    M_gov :: Vector{Float64}     # total amount of money at gov

    # debtt levels
    debt_tot :: Vector{Float64}   # total debt
    debt_cp :: Vector{Float64}    # cp debt
    debt_cp_allowed :: Vector{Float64}    # cp allowed debt
    debt_kp :: Vector{Float64}    # kp debt

    # Wage statistics
    w̄_avg :: Vector{Float64}     # average wage over time
    # w̄_std :: Vector{Float64}     # std of wage over time
    wʳ_avg :: Vector{Float64}   # average requested wage over time
    # wʳ_std :: Vector{Float64}   # std of requested wage over time
    wˢ_avg :: Vector{Float64}   # average of satisfying wage over time
    # wˢ_std :: Vector{Float64}   # std of satisfying wage over time
    wᴼ_max_mean :: Vector{Float64} # average of wᴼ_max

    Ī_avg :: Vector{Float64}     # average income over time
    Ī_std :: Vector{Float64}     # std of income over time
    B̄_avg :: Vector{Float64}     # average of budget over time
    B̄_std :: Vector{Float64}     # std of budget over time
    U :: Vector{Float64}         # unemployment over time
    Exp_UB :: Vector{Float64}    # total spending on UB
    AB :: Vector{Float64}        # average labor productivity over time
    l :: Vector{Float64}         # unfilled demand over time
    E_bar :: Vector{Float64}     # average competetiveness over time
    r :: Vector{Float64}         # interest rate over time
    Ls :: Vector{Float64}        # labor supply over time
    Ld :: Vector{Float64}        # labor demand over time

    # Changes in savings rate
    s̄_avg :: Vector{Float64}     # average savings rate over time
    s̄_std :: Vector{Float64}     # std of savings over time

    # Changes in labor demand
    ΔL̄_avg :: Vector{Float64}    # average desired labor change
    ΔL̄_std :: Vector{Float64}    # std desired labor change
    ΔL̄_cp_avg :: Vector{Float64} # average desired over time for cp
    ΔL̄_cp_std :: Vector{Float64} # std desired labor change for cp
    ΔL̄_kp_avg :: Vector{Float64} # ' ' for kp
    ΔL̄_kp_std :: Vector{Float64}

    # Investment levels
    EI_avg :: Vector{Float64}    # average expansion investment
    RS_avg :: Vector{Float64}    # average replacement investment

    # Productivity
    avg_π :: Vector{Float64}     # average productivity cp
    avg_A :: Vector{Float64}     # average A at kp
    avg_B :: Vector{Float64}     # average B at kp

    # Production
    avg_Q_bp :: Vector{Float64}  # average production of bp
    avg_Q_lp :: Vector{Float64}  # average production of lp
    avg_Q_kp :: Vector{Float64}  # average production of kp

end


function initialize_macro(
    T::Int
    )
    macro_struct = MacroEconomy(
        [],                     # GDP over time
        [],                     # income share of GDP
        [],                     # profit share of GDP of cp over time
        [],                     # profit share of GDP of kp over time
        [],                     # p̄: average price over time
        [],                     # CPI: inflation rate
        [],                     # aggregate consumption

        # Money amounts
        [],                     # M total
        [],                     # M hh
        [],                     # M cp
        [],                     # M kp
        [],                     # M gov

        # debtt levels
        [],                     # total debt
        [],                     # cp debt
        [],                     # cp allowed debt
        [],                     # kp debt

        [],                     # average of wage over time
        # [],                     # std of wage over time
        [],
        # [],
        [],
        # [],
        [],                     # wᴼ_max_mean: average of wᴼ_max

        [],                     # average income over time
        [],                     # std of income over time
        [],                     # average budget over time
        [],                     # std of budget over time
        [],                     
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],

        # Investments
        [],                     # EI
        [],                     # RS

        # Productivity
        [],
        [],
        [],

        # Production
        [],
        [],
        []
    )
    return macro_struct
end


"""
Updates macro stats after each time step
"""
function update_macro_timeseries(
    macro_struct::MacroEconomy, 
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    all_kp::Vector{Int},
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    E::Float64, 
    gov_struct::Government,
    global_param::GlobalParam,
    model::ABM
    )

    # Retrieve agent structs
    all_hh_str = map(hh_id -> model[hh_id], all_hh)
    all_cp_str = map(cp_id -> model[cp_id], all_cp)
    all_kp_str = map(kp_id -> model[kp_id], all_kp)
    all_bp_str = map(bp_id -> model[bp_id], all_bp)
    all_lp_str = map(lp_id -> model[lp_id], all_lp)

    # Compute GDP
    compute_GDP!(all_hh_str, all_cp_str, all_kp_str, macro_struct)

    # Compute CPI
    update_CPI!(all_cp_str, macro_struct)

    ΔL̄_avg = mean(map(p -> p.ΔLᵈ, vcat(all_cp_str, all_kp_str)))
    ΔL̄_std = std(map(p -> p.ΔLᵈ, vcat(all_cp_str, all_kp_str)))
    push!(macro_struct.ΔL̄_avg, ΔL̄_avg)
    push!(macro_struct.ΔL̄_std, ΔL̄_std)

    ΔL̄_cp_avg = mean(map(cp -> cp.ΔLᵈ, all_cp_str))
    ΔL̄_cp_std = std(map(cp -> cp.ΔLᵈ, all_cp_str))
    push!(macro_struct.ΔL̄_cp_avg, ΔL̄_cp_avg)
    push!(macro_struct.ΔL̄_cp_std, ΔL̄_cp_std)

    ΔL̄_kp_avg = mean(map(kp -> kp.ΔLᵈ, all_kp_str))
    ΔL̄_kp_std = std(map(kp -> kp.ΔLᵈ, all_kp_str))
    push!(macro_struct.ΔL̄_kp_avg, ΔL̄_kp_avg)
    push!(macro_struct.ΔL̄_kp_std, ΔL̄_kp_std)

    push!(macro_struct.U, E)

    push!(macro_struct.Exp_UB, gov_struct.curracc.Exp_UB[end])

    # Compute total amount in system
    compute_M!(all_hh_str, all_cp_str, all_kp_str, gov_struct, macro_struct)

    # Compute average savings rate
    s̄_avg = mean(map(hh -> hh.s, all_hh_str))
    s̄_std = std(map(hh -> hh.s, all_hh_str))
    # println(map(hh -> hh.s, all_hh_str))
    push!(macro_struct.s̄_avg, s̄_avg)
    push!(macro_struct.s̄_std, s̄_std)

    # Wage and income statistics
    update_wage_stats!(all_hh_str, all_cp_str, all_kp_str, macro_struct)

    Ī_avg = mean(map(hh -> hh.I[end], all_hh_str))
    Ī_std = std(map(hh -> hh.I[end], all_hh_str))
    push!(macro_struct.Ī_avg, Ī_avg)
    push!(macro_struct.Ī_std, Ī_std)

    D̄ = mean(map(cp -> cp.D[end], all_cp_str))
    # println(D̄)

    update_debt!(
        all_cp_str, 
        all_kp_str, 
        global_param.Λ, 
        macro_struct
    )

    # Investment
    EI_avg = mean(map(cp -> cp.EIᵈ, all_cp_str))
    push!(macro_struct.EI_avg, EI_avg)
    RS_avg = mean(map(cp -> cp.RSᵈ, all_cp_str))
    push!(macro_struct.RS_avg, RS_avg)

    # Productivity
    avg_π = mean(map(cp -> cp.π, all_cp_str))
    push!(macro_struct.avg_π, avg_π)
    avg_A = mean(map(kp -> kp.A[end], all_kp_str))
    push!(macro_struct.avg_A, avg_A)
    avg_B = mean(map(kp -> kp.B[end], all_kp_str))
    push!(macro_struct.avg_B, avg_B)

    # Production quantity
    avg_Q_bp = mean(map(bp -> bp.Q[end], all_bp_str))
    push!(macro_struct.avg_Q_bp, avg_Q_bp)
    # println(mean(map(cp -> cp.Qˢ, all_cp_str)))
    avg_Q_lp = mean(map(lp -> lp.Q[end], all_lp_str))
    push!(macro_struct.avg_Q_lp, avg_Q_lp)
    avg_Q_kp = mean(map(kp -> kp.Q[end], all_kp_str))
    push!(macro_struct.avg_Q_kp, avg_Q_kp)

end


"""
Computes GDP based on income of separate sectors, computes partial incomes of sectors
"""
function compute_GDP!(
    all_hh_str::Vector{Household},
    all_cp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    macro_struct::MacroEconomy
    )

    total_I = sum(map(hh -> hh.Iᵀ[end], all_hh_str))
    push!(macro_struct.GDP_I, total_I)

    total_Π_cp = sum(map(cp -> cp.Π[end], all_cp_str))

    push!(macro_struct.GDP_Π_cp, total_Π_cp)
    
    total_Π_kp = sum(map(kp -> kp.Π[end], all_kp_str))
    push!(macro_struct.GDP_Π_kp, total_Π_kp)

    # println("total I: ", total_I, ", Π cp: ", total_Π_cp, ", Π kp: ", total_Π_kp)

    GDP = total_I + total_Π_cp + total_Π_kp
    push!(macro_struct.GDP, GDP)
end

function update_labor_stats(macro_struct, labormarket_struct)

end

function compute_M!(
    all_hh_str::Vector{Household},
    all_cp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    gov_struct::Government,
    macro_struct::MacroEconomy
    )

    # Wealth of households
    M_hh = sum(map(hh -> hh.W[end], all_hh_str))
    # MA_avg = mean(map(hh -> hh.W[end], all_hh_str))
    # println("MA hh avg ", MA_avg)
    push!(macro_struct.M_hh, M_hh)

    # Liquid assets of cp
    M_cp = sum(map(cp -> cp.balance.NW - cp.balance.debt, all_cp_str))
    push!(macro_struct.M_cp, M_cp)

    # Liquid assets of kp
    M_kp = sum(map(kp -> kp.balance.NW - kp.balance.debt, all_kp_str))
    push!(macro_struct.M_kp, M_kp)

    # Money owned by government
    M_gov = gov_struct.MS
    push!(macro_struct.M_gov, M_gov)

    # Total amount of money stocks
    M_tot = M_hh + M_cp + M_kp + M_gov
    push!(macro_struct.M, M_tot)
end


function update_wage_stats!(
    all_hh_str::Vector{Household},
    all_cp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    macro_struct::MacroEconomy
    )

    w̄_avg = mean(map(p -> p.w̄[end], vcat(all_cp_str, all_kp_str)))
    # w̄_std = std(map(p -> p.w̄[end], vcat(all_cp_str, all_kp_str)))
    push!(macro_struct.w̄_avg, w̄_avg)
    # push!(macro_struct.w̄_std, w̄_std)

    wʳ_avg = mean(map(hh -> hh.wʳ, all_hh_str))
    # wʳ_std = std(map(hh -> hh.wʳ, all_hh_str))
    push!(macro_struct.wʳ_avg, wʳ_avg)
    # push!(macro_struct.wʳ_std, wʳ_std)

    wˢ_avg = mean(map(hh -> hh.wˢ, all_hh_str))
    # wˢ_std = std(map(hh -> hh.wˢ, all_hh_str))
    push!(macro_struct.wˢ_avg, wˢ_avg)
    # push!(macro_struct.wˢ_std, wˢ_std)

    wᴼ_max_mean = mean(map(p -> p.wᴼ_max, vcat(all_cp_str, all_kp_str)))
    push!(macro_struct.wᴼ_max_mean, wᴼ_max_mean)
end


"""
Updates metrics on aggregate debt levels
"""
function update_debt!(
    all_cp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    Λ::Float64,
    macro_struct::MacroEconomy
    )

    debt_cp = sum(map(cp -> cp.balance.debt, all_cp_str))
    push!(macro_struct.debt_cp, debt_cp)

    debt_kp = sum(map(kp -> kp.balance.debt, all_kp_str))
    push!(macro_struct.debt_kp, debt_kp)

    debt_tot = debt_cp + debt_kp
    push!(macro_struct.debt_tot, debt_tot)

    debt_cp_allowed = Λ * sum(map(cp -> cp.curracc.S, all_cp_str))
    push!(macro_struct.debt_cp_allowed, debt_cp_allowed)
end


function update_CPI!(
    all_cp_str::Vector{ConsumerGoodProducer},
    macro_struct::MacroEconomy
    )

    avg_p_1 = sum(map(cp -> cp.p[1] * cp.f[1], all_cp_str))
    avg_p_t = sum(map(cp -> cp.p[end] * cp.f[end], all_cp_str))

    push!(macro_struct.p̄, avg_p_t)
    # avg_w = mean(map(cp -> cp.w̄[end], all_cp_str))
    # avg_c = mean(map(cp -> cp.c[end], all_cp_str))
    # println(avg_p_1, " ", avg_p_t, " ", avg_w, " ", avg_c)
    # println(sum(map(cp -> 0.5*cp.f[end], all_cp_str)))

    # Compute average price, weighted by market share
    if length(macro_struct.CPI) == 0
        push!(macro_struct.CPI, 100)
    else
        # TODO dont compute this again every time
        CPI_t = 100 / avg_p_1 * avg_p_t
        push!(macro_struct.CPI, CPI_t)
    end
end