
mutable struct MacroEconomy
    GDP :: Vector{Float64}       # GDP over time
    GDP_I :: Vector{Float64}     # income share of GDP over time
    GDP_Π_cp :: Vector{Float64}  # profit share of GDP of cp over time
    GDP_Π_kp :: Vector{Float64}  # profit share of GDP of kp over time

    GDP_growth :: Vector{Float64}   # GDP growth rates over time

    p̄::Vector{Float64}           # average price of cp goods over time
    p̄_kp::Vector{Float64}        # average price of kp goods over time
    μ_bp::Vector{Float64}
    μ_lp::Vector{Float64}
    μ_kp::Vector{Float64}
    CPI :: Vector{Float64}       # price index consumer goods over time
    CPI_kp :: Vector{Float64}    # price index capital goods over time

    C :: Vector{Float64}         # aggregate consumption over time
    unsat_demand :: Vector{Float64} # average ratio of unsatisfied demand
    avg_N_goods :: Vector{Float64}

    # Division of money over sectors
    M :: Vector{Float64}         # total amount of money (should be stable)
    M_hh :: Vector{Float64}      # total amount of money at hh
    M_cp :: Vector{Float64}      # total amount of money at cp
    M_kp :: Vector{Float64}      # total amount of money at kp
    M_gov :: Vector{Float64}     # total amount of money at gov
    M_if :: Vector{Float64}      # total amount of money at if

    # debt levels
    debt_tot :: Vector{Float64}             # total debt
    debt_cp :: Vector{Float64}              # cp debt
    debt_cp_allowed :: Vector{Float64}      # cp allowed debt
    debt_kp :: Vector{Float64}              # kp debt
    debt_kp_allowed :: Vector{Float64}      # kp debt allowed
    debt_unpaid_kp :: Vector{Float64}       # kp debt that went unpaid after bankrupcy
    debt_unpaid_cp :: Vector{Float64}       # cp debt that went unpaid after bankrupcy

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
    EI_avg :: Vector{Float64}           # average expansion investment
    n_mach_EI_avg :: Vector{Float64}    # average amount of ordered machines for EI
    RS_avg :: Vector{Float64}           # average replacement investment
    n_mach_RS_avg :: Vector{Float64}    # average amounf of ordered machines for RS

    # Productivity
    avg_π :: Vector{Float64}     # average productivity cp
    avg_A :: Vector{Float64}     # average A at kp
    avg_B :: Vector{Float64}     # average B at kp

    # Production
    avg_Q_bp :: Vector{Float64}  # average production of bp
    avg_Q_lp :: Vector{Float64}  # average production of lp
    avg_Q_kp :: Vector{Float64}  # average production of kp

    # Bankrupties
    bankrupt_bp :: Vector{Float64}      # fraction of bp that went bankrupt
    bankrupt_lp :: Vector{Float64}      # fraction of lp that went bankrupt
    bankrupt_kp :: Vector{Float64}      # fraction of kp that went bankrupt

    cu :: Vector{Float64}               # average rate of capital utilization
    avg_n_machines_bp :: Vector{Float64}# average number of machines bp
    avg_n_machines_lp :: Vector{Float64}# average number of machines lp

    GINI_I::Vector{Float64}
    GINI_W::Vector{Float64}

end


function initialize_macro(
    T::Int
    )
    macro_struct = MacroEconomy(
        [],                     # GDP over time
        [],                     # income share of GDP
        [],                     # profit share of GDP of cp over time
        [],                     # profit share of GDP of kp over time

        [],                     # GDP growth rates over time

        [],                     # p̄: average price over time
        [],                     # p̄_kp: average price of kp goods over time
        [],
        [],
        [],
        [],                     # CPI: inflation rate
        [],                     # CPI_kp: price index capital goods

        [],                     # aggregate consumption
        [],                     # unsat_dem: ratio of unsatisfied demand
        [],                     # Inventory ratio

        # Money amounts
        [],                     # M total
        [],                     # M hh
        [],                     # M cp
        [],                     # M kp
        [],                     # M gov
        [],                     # M if

        # debt levels
        [],                     # total debt
        [],                     # cp debt
        [],                     # cp allowed debt
        [],                     # kp debt
        [],                     # kp allowed debt
        [],                     # kp unpaid debt
        [],                     # cp unpaid debt

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
        [],                     # n mach EI
        [],                     # RS
        [],                     # n mach RS

        # Productivity
        [],
        [],
        [],

        # Production
        [],
        [],
        [],

        # Bankrupties
        [],
        [],
        [],

        [],                     # cu: capital utilization rate
        [],                     # avg_n_machines_bp
        [],                     # avg_n_machines_lp

        # Inequality
        [],                     # GINI_I
        [],                     # GINI_W
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
    bankrupt_bp::Vector{Int},
    bankrupt_lp::Vector{Int},
    bankrupt_kp::Vector{Int},
    E::Float64, 
    gov_struct::Government,
    indexfund_struct::IndexFund,
    global_param::GlobalParam,
    model::ABM
    )

    # Retrieve agent structs
    all_hh_str = map(hh_id -> model[hh_id], all_hh)
    all_cp_str = map(cp_id -> model[cp_id], all_cp)
    all_kp_str = map(kp_id -> model[kp_id], all_kp)
    all_bp_str = map(bp_id -> model[bp_id], all_bp)
    all_lp_str = map(lp_id -> model[lp_id], all_lp)


    avg_μ = mean(map(p -> p.μ[end], vcat(all_cp_str)))
    # println("avg μ: $avg_μ")

    # Compute GDP
    compute_GDP!(all_hh_str, all_cp_str, all_kp_str, macro_struct)

    # Compute CPI
    compute_price_data!(
        all_cp_str,
        all_bp_str,
        all_lp_str, 
        all_kp_str, 
        macro_struct
    )

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

    # Consumption 
    C = sum(map(hh->hh.C[end], all_hh_str))
    push!(macro_struct.C, C)

    push!(macro_struct.U, E)

    push!(macro_struct.Exp_UB, gov_struct.curracc.Exp_UB[end])

    # Compute total amount in system
    compute_M!(all_hh_str, all_cp_str, all_kp_str, gov_struct, 
               indexfund_struct, macro_struct)

    # Compute average savings rate
    s̄_avg = mean(map(hh -> hh.s, all_hh_str))
    s̄_std = std(map(hh -> hh.s, all_hh_str))
    push!(macro_struct.s̄_avg, s̄_avg)
    push!(macro_struct.s̄_std, s̄_std)

    # Wage and income statistics
    update_wage_stats!(all_hh_str, all_cp_str, all_kp_str, macro_struct)

    Ī_avg = mean(map(hh -> hh.I[end], all_hh_str))
    Ī_std = std(map(hh -> hh.I[end], all_hh_str))
    push!(macro_struct.Ī_avg, Ī_avg)
    push!(macro_struct.Ī_std, Ī_std)

    # D̄ = mean(map(cp -> cp.D[end], all_cp_str))
    # println(D̄)

    update_debt!(
        all_cp_str, 
        all_kp_str,
        bankrupt_bp,
        bankrupt_lp, 
        bankrupt_kp,
        global_param.Λ, 
        macro_struct,
        model
    )

    # Investment
    EI_avg = mean(map(cp -> cp.EIᵈ, all_cp_str))
    push!(macro_struct.EI_avg, EI_avg)
    n_mach_EI_avg = mean(cp -> cp.n_mach_ordered_EI, all_cp_str)
    push!(macro_struct.n_mach_EI_avg, n_mach_EI_avg)
    RS_avg = mean(map(cp -> cp.RSᵈ, all_cp_str))
    push!(macro_struct.RS_avg, RS_avg)
    n_mach_RS_avg = mean(cp -> cp.n_mach_ordered_RS, all_cp_str)
    push!(macro_struct.n_mach_RS_avg, n_mach_RS_avg)

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

    compute_bankrupties(
        all_bp, 
        all_lp, 
        all_kp, 
        bankrupt_bp, 
        bankrupt_lp, 
        bankrupt_kp,
        macro_struct
    )

    compute_unsatisfied_demand(
        all_hh_str,
        macro_struct
    )

    avg_N_goods = mean(map(cp -> cp.N_goods, all_cp_str))
    push!(macro_struct.avg_N_goods, avg_N_goods)

    # Mean rate of capital utilization
    cu = mean(map(cp -> cp.cu, all_cp_str))
    push!(macro_struct.cu, cu)

    avg_n_machines_bp = mean(map(bp->bp.n_machines, all_bp_str))
    push!(macro_struct.avg_n_machines_bp, avg_n_machines_bp)

    avg_n_machines_lp = mean(map(lp->lp.n_machines, all_lp_str))
    push!(macro_struct.avg_n_machines_lp, avg_n_machines_lp)

    # println("mean D: ", mean(map(cp->cp.D[end], all_cp_str)),
    #         " mean Dᵉ: ", mean(map(cp->cp.Dᵉ, all_cp_str)),
    #         " mean Qˢ: ", mean(map(cp->cp.Qˢ, all_cp_str)))

    # println("sum D: ", sum(map(cp->cp.D[end], all_cp_str)), 
    #         " sum N_goods: ", sum(map(hh->hh.N_goods, all_hh_str)),
    #         " sum Q: ", sum(map(cp->cp.Q[end], all_cp_str)))

    # Compute GINI coefficients
    compute_GINI(all_hh_str, macro_struct)

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

    if length(GDP) > 1
        GDP_growth = (GDP[end] - GDP[end-1]) / GDP[end-1]
    else
        GDP_growth = 0
    end
    push!(macro_struct.GDP_growth, GDP_growth)
end

function update_labor_stats(macro_struct, labormarket_struct)

end


"""
Computes the ratios of bankrupt bp, lp and kp.
"""
function compute_bankrupties(
    all_bp,
    all_lp,
    all_kp,
    bankrupt_bp,
    bankrupt_lp,
    bankrupt_kp,
    macro_struct
    )

    bankrupt_bp = length(bankrupt_bp) / length(all_bp)
    push!(macro_struct.bankrupt_bp, bankrupt_bp)

    bankrupt_lp = length(bankrupt_lp) / length(all_lp)
    push!(macro_struct.bankrupt_lp, bankrupt_lp)

    bankrupt_kp = length(bankrupt_kp) / length(all_kp)
    push!(macro_struct.bankrupt_kp, bankrupt_kp)

end


function compute_M!(
    all_hh_str::Vector{Household},
    all_cp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    gov_struct::Government,
    indexfund_struct::IndexFund,
    macro_struct::MacroEconomy
    )

    # Wealth of households
    M_hh = sum(map(hh -> hh.W[end], all_hh_str))
    # MA_avg = mean(map(hh -> hh.W[end], all_hh_str))
    # println("MA hh avg ", MA_avg)
    push!(macro_struct.M_hh, M_hh)

    # Liquid assets of cp
    M_cp = sum(map(cp -> cp.balance.NW, all_cp_str))
    push!(macro_struct.M_cp, M_cp)

    # Liquid assets of kp
    M_kp = sum(map(kp -> kp.balance.NW, all_kp_str))
    push!(macro_struct.M_kp, M_kp)

    # Money owned by government
    M_gov = gov_struct.MS
    push!(macro_struct.M_gov, M_gov)

    # Money in investment fund
    push!(macro_struct.M_if, indexfund_struct.Assets)

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
    bankrupt_bp::Vector{Int},
    bankrupt_lp::Vector{Int},
    bankrupt_kp::Vector{Int},
    Λ::Float64,
    macro_struct::MacroEconomy,
    model::ABM
    )

    debt_cp = sum(map(cp -> cp.balance.debt, all_cp_str))
    push!(macro_struct.debt_cp, debt_cp)

    debt_kp = sum(map(kp -> kp.balance.debt, all_kp_str))
    push!(macro_struct.debt_kp, debt_kp)

    debt_tot = debt_cp + debt_kp
    push!(macro_struct.debt_tot, debt_tot)

    debt_cp_allowed = Λ * sum(map(cp -> cp.curracc.S, all_cp_str))
    push!(macro_struct.debt_cp_allowed, debt_cp_allowed)

    debt_kp_allowed = Λ * sum(map(kp -> kp.curracc.S, all_kp_str))
    push!(macro_struct.debt_kp_allowed, debt_kp_allowed)

    debt_unpaid_cp = sum(map(cp_id -> model[cp_id].balance.debt, vcat(bankrupt_bp, bankrupt_lp)))
    push!(macro_struct.debt_unpaid_cp, debt_unpaid_cp)

    debt_unpaid_kp = sum(map(kp_id -> model[kp_id].balance.debt, bankrupt_kp))
    push!(macro_struct.debt_unpaid_kp, debt_unpaid_kp)
end


function compute_price_data!(
    all_cp_str::Vector{ConsumerGoodProducer},
    all_bp_str::Vector{ConsumerGoodProducer},
    all_lp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    macro_struct::MacroEconomy
    )

    # Compute average price, weighted by market share
    avg_p_t = sum(map(cp -> cp.p[end] * cp.f[end], all_cp_str))
    push!(macro_struct.p̄, avg_p_t)
    if length(macro_struct.CPI) == 0
        push!(macro_struct.CPI, 100)
    else
        # TODO dont compute this again every time
        CPI_t = 100 / macro_struct.p̄[1] * avg_p_t
        push!(macro_struct.CPI, CPI_t)
    end

    # Compute average price of capital goods
    avg_p_kp_t = mean(map(kp -> kp.p[end], all_kp_str))
    push!(macro_struct.p̄_kp, avg_p_kp_t)
    if length(macro_struct.CPI_kp) == 0
        push!(macro_struct.CPI_kp, 100)
    else
        # TODO dont compute this again every time
        CPI_kp_t = 100 / macro_struct.p̄_kp[1] * avg_p_kp_t
        push!(macro_struct.CPI_kp, CPI_kp_t)
    end

    # Update markup rates
    μ_bp = mean(map(bp -> bp.μ[end], all_bp_str))
    push!(macro_struct.μ_bp, μ_bp)

    μ_lp = mean(map(lp -> lp.μ[end], all_lp_str))
    push!(macro_struct.μ_lp, μ_lp)

    μ_kp = mean(map(kp -> kp.μ[end], all_kp_str))
    push!(macro_struct.μ_kp, μ_kp)

end


function compute_unsatisfied_demand(
    all_hh_str::Vector{Household},
    macro_struct::MacroEconomy
    )

    avg_unsat_dem = Vector{Float64}()

    for hh in all_hh_str
        # println(hh.unsat_dem)
        if length(hh.unsat_dem) > 0 
            unsat_dem = mean(map(ud->ud[2], hh.unsat_dem))
        else
            unsat_dem = 0.0
        end
        push!(avg_unsat_dem, unsat_dem)
    end

    avg_unsat_dem = mean(avg_unsat_dem)
    push!(macro_struct.unsat_demand, avg_unsat_dem)
end


"""
Computes the GINI coefficient for wealth and income
"""
function compute_GINI(
    all_hh_str::Vector{Household},
    macro_struct::MacroEconomy
    )

    # Compute GINI for income
    all_I = map(hh -> hh.Iᵀ[end], all_hh_str)
    GINI_I = sum(map(hh -> sum(broadcast(abs, hh.Iᵀ[end] .- all_I)), all_hh_str))
    GINI_I = GINI_I / (2 * length(all_hh_str)^2 * macro_struct.Ī_avg[end])
    push!(macro_struct.GINI_I, GINI_I)

    # Compute GINI for wealth
    all_W = map(hh -> hh.W[end], all_hh_str)
    GINI_W = sum(map(hh -> sum(broadcast(abs, hh.W[end] .- all_W)), all_hh_str))
    GINI_W = GINI_W / (2 * length(all_hh_str)^2 * macro_struct.M_hh[end] / length(all_hh_str))
    push!(macro_struct.GINI_W, GINI_W)

    # println("GINI income: $GINI_I, GINI wealth: $GINI_W")

end