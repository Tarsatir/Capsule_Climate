
@Base.kwdef mutable struct MacroEconomy
    T::Int=T                                                # number of timesteps

    # GDP
    GDP::Vector{Float64} = zeros(Float64, T)                # GDP over time
    GDP_I::Vector{Float64} = zeros(Float64, T)              # income share of GDP over time
    GDP_Π_cp::Vector{Float64} = zeros(Float64, T)           # profit share of GDP of cp over time
    GDP_Π_kp::Vector{Float64} = zeros(Float64, T)           # profit share of GDP of kp over time
    GDP_growth::Vector{Float64} = zeros(Float64, T)         # quarterly GDP growth rates

    total_C::Vector{Float64} = zeros(Float64, T)            # total spending on consumption
    total_C_actual::Vector{Float64} = zeros(Float64, T)     # total actual spending on consumption
    total_I::Vector{Float64} = zeros(Float64, T)            # total actual spending on investments
    total_w::Vector{Float64} = zeros(Float64, T)            # total actual spending on wages

    LIS::Vector{Float64} = zeros(Float64, T)                # labor income share

    returns_investments::Vector{Float64} = zeros(Float64, T)# returns to investments

    p̄::Vector{Float64} = zeros(Float64, T)                  # average price of cp goods over time
    p̄_kp::Vector{Float64} = zeros(Float64, T)               # average price of kp goods over time
    μ_cp::Vector{Float64} = zeros(Float64, T)               # average μ for cp
    μ_kp::Vector{Float64} = zeros(Float64, T)               # average μ for kp
    CPI::Vector{Float64} = zeros(Float64, T)                # price index consumer goods over time
    CPI_kp::Vector{Float64} = zeros(Float64, T)             # price index capital goods over time

    unsat_demand::Vector{Float64} = zeros(Float64, T)       # ratio of unsatisfied demand
    unsat_invest::Vector{Float64} = zeros(Float64, T)       # ratio of unsatisfied investments
    unsat_L_demand::Vector{Float64} = zeros(Float64, T)     # ratio of unsatisfied labor demand
    unspend_C::Vector{Float64} = zeros(Float64, T)          # average ratio of unspend C
    N_goods::Vector{Float64} = zeros(Float64, T)            # total amount of inventories
    avg_N_goods::Vector{Float64} = zeros(Float64, T)        # average number of goods in cp firm inventory

    # Division of money over sectors
    M::Vector{Float64} = zeros(Float64, T)                  # total amount of money (should be stable)
    M_hh::Vector{Float64} = zeros(Float64, T)               # total amount of money at hh
    M_cp::Vector{Float64} = zeros(Float64, T)               # total amount of money at cp
    M_kp::Vector{Float64} = zeros(Float64, T)               # total amount of money at kp
    M_ep::Vector{Float64} = zeros(Float64, T)               # total amount of money at ep
    M_gov::Vector{Float64} = zeros(Float64, T)              # total amount of money at gov
    M_if::Vector{Float64}  = zeros(Float64, T)              # total amount of money at if

    # debt levels
    debt_tot::Vector{Float64} = zeros(Float64, T)           # total debt
    debt_cp::Vector{Float64} = zeros(Float64, T)            # cp debt
    debt_cp_allowed::Vector{Float64} = zeros(Float64, T)    # cp allowed debt
    debt_kp::Vector{Float64} = zeros(Float64, T)            # kp debt
    debt_kp_allowed::Vector{Float64} = zeros(Float64, T)    # kp debt allowed
    debt_unpaid_kp::Vector{Float64} = zeros(Float64, T)     # kp debt that went unpaid after bankrupcy
    debt_unpaid_cp::Vector{Float64} = zeros(Float64, T)     # cp debt that went unpaid after bankrupcy

    # Wage statistics
    w̄_avg::Vector{Float64} = ones(Float64, T)               # average wage over time
    wʳ_avg::Vector{Float64} = ones(Float64, T)              # average requested wage over time
    wˢ_avg::Vector{Float64} = ones(Float64, T)              # average of satisfying wage over time
    wᴼ_max_mean::Vector{Float64} = ones(Float64, T)         # average of wᴼ_max

    Ī_avg::Vector{Float64} = zeros(Float64, T)              # average income over time
    I_labor_avg::Vector{Float64} = zeros(Float64, T)        # average labor income over time
    I_capital_avg::Vector{Float64} = zeros(Float64, T)      # average capital income over time
    I_UB_avg::Vector{Float64} = zeros(Float64, T)           # average UB income over time
    I_socben_avg::Vector{Float64} = zeros(Float64, T)       # average social benefits income over time

    B̄_avg::Vector{Float64} = zeros(Float64, T)              # average of budget over time
    B̄_std::Vector{Float64} = zeros(Float64, T)              # std of budget over time

    U::Vector{Float64} = zeros(Float64, T)                  # unemployment over time
    dU::Vector{Float64} = zeros(Float64, T)                 # quarterly change in unemployment rate
    switch_rate::Vector{Float64} = zeros(Float64, T)        # rate of switching employers
    Exp_UB::Vector{Float64} = zeros(Float64, T)             # total spending on UB
    AB::Vector{Float64} = zeros(Float64, T)                 # average labor productivity over time
    l::Vector{Float64} = zeros(Float64, T)                  # unfilled demand over time
    E_bar::Vector{Float64} = zeros(Float64, T)              # average competetiveness over time
    r::Vector{Float64} = zeros(Float64, T)                  # interest rate over time
    Ls::Vector{Float64} = zeros(Float64, T)                 # labor supply over time
    Ld::Vector{Float64} = zeros(Float64, T)                 # labor demand over time

    # Changes in savings rate
    s̄_emp::Vector{Float64} = zeros(Float64, T)              # average savings rate of employed households over time
    s̄_unemp::Vector{Float64} = zeros(Float64, T)            # average savings rate of unemployed households over time

    # Changes in labor demand
    ΔL̄_avg::Vector{Float64} = zeros(Float64, T)             # average desired labor change
    ΔL̄_std::Vector{Float64} = zeros(Float64, T)             # std desired labor change
    ΔL̄_cp_avg::Vector{Float64} = zeros(Float64, T)          # average desired over time for cp
    ΔL̄_cp_std::Vector{Float64} = zeros(Float64, T)          # std desired labor change for cp
    ΔL̄_kp_avg::Vector{Float64} = zeros(Float64, T)          # average desired over time for kp
    ΔL̄_kp_std::Vector{Float64} = zeros(Float64, T)          # std desired labor change for kp

    # Investment levels
    RD_total::Vector{Float64} = zeros(Float64, T)           # total R&D investments
    EI_avg::Vector{Float64}  = zeros(Float64, T)            # average expansion investment
    n_mach_EI_avg::Vector{Float64} = zeros(Float64, T)      # average amount of ordered machines for EI
    RS_avg::Vector{Float64} = zeros(Float64, T)             # average replacement investment
    n_mach_RS_avg::Vector{Float64} = zeros(Float64, T)      # average amounf of ordered machines for RS

    # Productivity
    avg_π_LP::Vector{Float64} = zeros(Float64, T)           # average labor productivity cp
    avg_π_EE::Vector{Float64} = zeros(Float64, T)           # average productivity per energy unit cp
    avg_π_EF::Vector{Float64} = zeros(Float64, T)           # average energy friendliness 

    avg_A_LP::Vector{Float64} = zeros(Float64, T)           # average A_LP at kp
    avg_A_EE::Vector{Float64} = zeros(Float64, T)           # average A_EE at kp
    avg_A_EF::Vector{Float64} = zeros(Float64, T)           # average A_EF at kp

    avg_B_LP::Vector{Float64} = zeros(Float64, T)           # average B_LP at kp
    avg_B_EE::Vector{Float64} = zeros(Float64, T)           # average B_EE at kp
    avg_B_EF::Vector{Float64} = zeros(Float64, T)           # average B_EF at kp

    # Production
    total_Q_cp::Vector{Float64} = zeros(Float64, T)         # total units produced by cp
    total_Q_kp::Vector{Float64} = zeros(Float64, T)         # total units produced by kp
    total_Q_growth::Vector{Float64} = zeros(Float64, T)     # growth in total units produced

    avg_Q_cp::Vector{Float64} = zeros(Float64, T)           # average production of cp
    avg_Qˢ_cp::Vector{Float64} = zeros(Float64, T)          # average desired ST production of cp
    avg_Qᵉ_cp::Vector{Float64} = zeros(Float64, T)          # average desired LT production of cp
    avg_Q_kp::Vector{Float64} = zeros(Float64, T)           # average production of kp
    avg_D_cp::Vector{Float64} = zeros(Float64, T)           # average demand of cp
    avg_Dᵁ_cp::Vector{Float64} = zeros(Float64, T)          # average unsatisfied demand of cp
    avg_Dᵉ_cp::Vector{Float64} = zeros(Float64, T)          # average expected demand of cp

    # Bankrupties
    bankrupt_cp::Vector{Float64} = zeros(Float64, T)        # fraction of cp that went bankrupt
    bankrupt_kp::Vector{Float64} = zeros(Float64, T)        # fraction of kp that went bankrupt

    cu::Vector{Float64} = zeros(Float64, T)                 # average rate of capital utilization
    avg_n_machines_cp::Vector{Float64} = zeros(Float64, T)  # average number of machines cp

    GINI_I::Vector{Float64} = zeros(Float64, T)             # Gini coefficient for income
    GINI_W::Vector{Float64} = zeros(Float64, T)             # Gini coefficient for wealth

    I_min::Vector{Float64} = zeros(Float64, T)              # Minimum income in the economy
    I_20::Vector{Float64} = zeros(Float64, T)               # Threshold of the lower 20% of income
    I_80::Vector{Float64} = zeros(Float64, T)               # Threshold of the lower 80% of income
    I_max::Vector{Float64} = zeros(Float64, T)              # Maximum income in the economy

    W_min::Vector{Float64} = zeros(Float64, T)              # Minimum wealth level in the economy
    W_20::Vector{Float64} = zeros(Float64, T)               # Threshold of the lower 20% of wealth
    W_80::Vector{Float64} = zeros(Float64, T)               # Threshold of the lower 80% of wealth
    W_max::Vector{Float64} = zeros(Float64, T)              # Maximum wealth level in the economy

    α_W_quantiles::Matrix{Float64} = zeros(Float64, 5, T)   # Matrix containing average α for different wealth quintiles
end


"""
Updates macro stats after each time step
"""
function update_macro_timeseries(
    # macroeconomy::MacroEconomy,
    t::Int, 
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    all_kp::Vector{Int},
    all_p::Vector{Int},
    ep,
    bankrupt_cp::Vector{Int},
    bankrupt_kp::Vector{Int},
    labormarket, 
    government::Government,
    indexfund::IndexFund,
    globalparam::GlobalParam,
    model::ABM,
    to
    )

    # Update CPI
    # compute_price_data!(all_cp, all_kp, t, model.macroeconomy, model)
    compute_price_data!(model)

    # Update total GDP, per sector and GDP growth
    # compute_GDP!(all_hh, all_cp, all_kp, t, model)
    compute_GDP!(model)

    # Update spending of different sectors
    # compute_spending!(all_hh, all_cp, all_kp, all_p, t, model)
    compute_spending!(model)

    # Compute the labor share of income
    model.macroeconomy.LIS[t] = model.macroeconomy.total_w[t] / model.macroeconomy.GDP[t]

    # Update returns to investments
    model.macroeconomy.returns_investments[t] = indexfund.returns_investments

    # Update average labor demand
    model.macroeconomy.ΔL̄_avg[t] = mean(p_id -> model[p_id].ΔLᵈ, all_p)
    model.macroeconomy.ΔL̄_cp_avg[t] = mean(cp_id -> model[cp_id].ΔLᵈ, all_cp)
    model.macroeconomy.ΔL̄_kp_avg[t] = mean(kp_id -> model[kp_id].ΔLᵈ, all_kp)

    # Update unemployment rate and unemployment benefits expenditure
    model.macroeconomy.U[t] = labormarket.E
    model.macroeconomy.switch_rate[t] = labormarket.switch_rate
    model.macroeconomy.Exp_UB[t] = government.curracc.Exp_UB[t]

    # Compute change in unemployment rate
    if t > 3 
        model.macroeconomy.dU[t] = 100 * (model.macroeconomy.U[t] - model.macroeconomy.U[t-3]) / model.macroeconomy.U[t-3]
    end

    # Consumption 
    # model.macroeconomy.C[t] = sum(hh_id->model[hh_id].C[end], all_hh)

    # Compute total amount in system
    compute_M!(all_hh, all_cp, all_kp, ep, government, indexfund, t, model)

    # Compute average savings rates
    compute_savings_macro!(all_hh, t, model)

    # Wage and income statistics
    update_wage_stats!(all_hh, all_p, t, model)
    update_income_stats!(all_hh, t, model)

    update_debt!(all_cp, all_kp, bankrupt_cp, bankrupt_kp, globalparam.Λ, t, model)

    # Investment
    model.macroeconomy.RD_total[t] = sum(kp_id -> model[kp_id].RD, all_kp) + ep.RDₑ[t]
    model.macroeconomy.EI_avg[t] = mean(cp_id -> model[cp_id].EIᵈ, all_cp)
    model.macroeconomy.n_mach_EI_avg[t] = mean(cp_id -> model[cp_id].n_mach_ordered_EI, all_cp)
    model.macroeconomy.RS_avg[t] = mean(cp_id -> model[cp_id].RSᵈ, all_cp)
    model.macroeconomy.n_mach_RS_avg[t] = mean(cp_id -> model[cp_id].n_mach_ordered_RS, all_cp)

    # Productivity
    model.macroeconomy.avg_π_LP[t] = mean(cp_id -> model[cp_id].π_LP, all_cp)
    model.macroeconomy.avg_π_EE[t] = mean(cp_id -> model[cp_id].π_EE, all_cp)
    model.macroeconomy.avg_π_EF[t] = mean(cp_id -> model[cp_id].π_EF, all_cp)

    model.macroeconomy.avg_A_LP[t] = mean(kp_id -> model[kp_id].A_LP[end], all_kp)
    model.macroeconomy.avg_A_EE[t] = mean(kp_id -> model[kp_id].A_EE[end], all_kp)
    model.macroeconomy.avg_A_EF[t] = mean(kp_id -> model[kp_id].A_EF[end], all_kp)

    model.macroeconomy.avg_B_LP[t] = mean(kp_id -> model[kp_id].B_LP[end], all_kp)
    model.macroeconomy.avg_B_EE[t] = mean(kp_id -> model[kp_id].B_EE[end], all_kp)
    model.macroeconomy.avg_B_EF[t] = mean(kp_id -> model[kp_id].B_EF[end], all_kp)

    # Production quantity
    model.macroeconomy.total_Q_cp[t] = sum(cp_id -> model[cp_id].Q[end], all_cp)
    model.macroeconomy.total_Q_kp[t] = sum(kp_id -> model[kp_id].Q[end], all_kp)
    if t > 3
        total_Q_t = model.macroeconomy.total_Q_cp[t] + model.macroeconomy.total_Q_kp[t]
        total_Q_t3 = model.macroeconomy.total_Q_cp[t-3] + model.macroeconomy.total_Q_kp[t-3]
        model.macroeconomy.total_Q_growth[t] = (total_Q_t - total_Q_t3) / total_Q_t3
    end

    model.macroeconomy.avg_Q_cp[t] = mean(cp_id -> model[cp_id].Q[end], all_cp)
    model.macroeconomy.avg_Qˢ_cp[t] = mean(cp_id -> model[cp_id].Qˢ, all_cp)
    model.macroeconomy.avg_Qᵉ_cp[t] = mean(cp_id -> model[cp_id].Qᵉ, all_cp)
    model.macroeconomy.avg_Q_kp[t] = mean(kp_id -> model[kp_id].Q[end], all_kp)
    model.macroeconomy.avg_D_cp[t] = mean(cp_id -> model[cp_id].D[end], all_cp)
    model.macroeconomy.avg_Dᵁ_cp[t] = mean(cp_id -> model[cp_id].Dᵁ[end], all_cp)
    model.macroeconomy.avg_Dᵉ_cp[t] = mean(cp_id -> model[cp_id].Dᵉ, all_cp)

    compute_bankrupties!(all_cp, all_kp, bankrupt_cp, bankrupt_kp, t, model)

    compute_unsatisfied_demand(all_cp, all_kp, all_hh, labormarket, 
                               globalparam.freq_per_machine, t, model)

    model.macroeconomy.N_goods[t] = sum(cp_id -> model[cp_id].N_goods, all_cp)
    model.macroeconomy.avg_N_goods[t] = mean(cp_id -> model[cp_id].N_goods, all_cp)

    # Mean rate of capital utilization
    model.macroeconomy.cu[t] = mean(cp_id -> model[cp_id].n_machines > 0 ? model[cp_id].cu : 0.5, all_cp)

    # Average number of machines
    model.macroeconomy.avg_n_machines_cp[t] = mean(cp_id -> model[cp_id].n_machines, all_cp)

    # Compute GINI coefficients
    @timeit to "GINI" compute_GINI(all_hh, t, model)

    compute_I_W_thresholds(all_hh, t, model)

    compute_α_W_quantiles(all_hh, t, model)

end


"""
Computes GDP based on income of separate sectors, computes partial incomes of sectors
"""
function compute_GDP!(
    # all_hh::Vector{Int},
    # all_cp::Vector{Int},
    # all_kp::Vector{Int},
    # model.macroeconomy::MacroEconomy,
    # t::Int,
    model::ABM
    )

    # Household income
    total_I = sum(hh_id -> model[hh_id].total_I[end], model.all_hh)
    model.macroeconomy.GDP_I[model.t] = total_I

    # cp profits
    total_Π_cp = sum(cp_id -> model[cp_id].Π[end], model.all_cp)
    model.macroeconomy.GDP_Π_cp[model.t] = total_Π_cp
    
    # kp profits
    total_Π_kp = sum(kp_id -> model[kp_id].Π[end], model.all_kp)
    model.macroeconomy.GDP_Π_kp[model.t] = total_Π_kp

    # total GDP
    model.macroeconomy.GDP[model.t] = total_I + total_Π_cp + total_Π_kp

    # quarterly real GDP growth
    if model.t > 3
        real_GDP_t = model.macroeconomy.GDP[model.t] / model.macroeconomy.p̄[model.t]
        real_GDP_t3 = model.macroeconomy.GDP[model.t-3] / model.macroeconomy.p̄[model.t-3]
        model.macroeconomy.GDP_growth[model.t] = 100 * (real_GDP_t - real_GDP_t3) / real_GDP_t3
    end
end


function compute_spending!(
    # all_hh::Vector{Int64}, 
    # all_cp::Vector{Int64}, 
    # all_kp::Vector{Int64}, 
    # all_p::Vector{Int64},
    # model.macroeconomy::MacroEconomy, 
    # t::Int, 
    model::ABM
    )

    # Compute planned and actual consumption
    model.macroeconomy.total_C[model.t] = sum(hh_id -> model[hh_id].C, model.all_hh)
    model.macroeconomy.total_C_actual[model.t] = sum(hh_id -> model[hh_id].C_actual, model.all_hh)

    # Compute total spending on investments
    model.macroeconomy.total_I[model.t] = sum(cp_id -> model[cp_id].curracc.TCI, model.all_cp)
    
    # Compute total spending on wages
    model.macroeconomy.total_w[model.t] = sum(p_id -> model[p_id].curracc.TCL, model.all_p)
end


"""
Computes the ratios of bankrupt bp, lp and kp.
"""
function compute_bankrupties!(
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    bankrupt_cp::Vector{Int64},
    bankrupt_kp::Vector{Int64},
    # model.macroeconomy::MacroEconomy,
    t::Int64,
    model::ABM
    )

    model.macroeconomy.bankrupt_cp[t] = length(bankrupt_cp) / length(all_cp)

    model.macroeconomy.bankrupt_kp[t] = length(bankrupt_kp) / length(all_kp)

end


"""
Computes the money supply of households, producers, government 
    and the indexfund
"""
function compute_M!(
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    ep,
    government::Government,
    indexfund::IndexFund,
    # model.macroeconomy::MacroEconomy,
    t::Int,
    model::ABM
    )

    # Wealth of households
    model.macroeconomy.M_hh[t] = sum(hh_id -> model[hh_id].W, all_hh)

    # Liquid assets of cp_id
    model.macroeconomy.M_cp[t] = sum(cp_id -> model[cp_id].balance.NW, all_cp)

    # Liquid assets of kp
    model.macroeconomy.M_kp[t] = sum(kp_id -> model[kp_id].balance.NW, all_kp)

    # Liquid assets of ep
    model.macroeconomy.M_ep[t] = ep.NWₑ[t]

    # Money owned by government
    model.macroeconomy.M_gov[t] = government.MS

    # Money in investment fund
    model.macroeconomy.M_if[t] = indexfund.Assets

    # Total amount of money stocks
    model.macroeconomy.M[t] = (model.macroeconomy.M_hh[t] + model.macroeconomy.M_cp[t] + model.macroeconomy.M_kp[t] +
                         model.macroeconomy.M_ep[t] + model.macroeconomy.M_gov[t] + model.macroeconomy.M_if[t])
end


"""
Computes average wage statistics
"""
function update_wage_stats!(
    all_hh::Vector{Int},
    all_p::Vector{Int},
    # model.macroeconomy::MacroEconomy,
    t::Int,
    model::ABM
    )

    model.macroeconomy.w̄_avg[t] = mean(p_id -> model[p_id].w̄[end], all_p)

    model.macroeconomy.wʳ_avg[t] = mean(hh_id -> model[hh_id].wʳ, all_hh)

    model.macroeconomy.wˢ_avg[t] = mean(hh_id -> model[hh_id].wˢ, all_hh)

    model.macroeconomy.wᴼ_max_mean[t] = mean(p_id -> model[p_id].wᴼ_max, all_p)
end

function update_income_stats!(
    all_hh::Vector{Int}, 
    # model.macroeconomy::MacroEconomy, 
    t::Int, 
    model::ABM
    )

    model.macroeconomy.I_labor_avg[t] = mean(hh_id -> model[hh_id].labor_I, all_hh)

    model.macroeconomy.I_capital_avg[t] = mean(hh_id -> model[hh_id].capital_I, all_hh)

    model.macroeconomy.I_UB_avg[t] = mean(hh_id -> model[hh_id].UB_I, all_hh)

    model.macroeconomy.I_socben_avg[t] = mean(hh_id -> model[hh_id].socben_I, all_hh)

    model.macroeconomy.Ī_avg[t] = (model.macroeconomy.I_labor_avg[t] + model.macroeconomy.I_capital_avg[t] 
                             + model.macroeconomy.I_UB_avg[t] + model.macroeconomy.I_socben_avg[t])
    
    # all_I = map(hh_id -> model[hh_id].total_I, all_hh)
    # println(percentile(all_I, 0.1))
end


"""
Updates metrics on aggregate debt levels
"""
function update_debt!(
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    bankrupt_cp::Vector{Int},
    bankrupt_kp::Vector{Int},
    Λ::Float64,
    # macroeconomy::MacroEconomy,
    t::Int,
    model::ABM
    )

    model.macroeconomy.debt_cp[t] = sum(cp_id -> model[cp_id].balance.debt, all_cp)

    model.macroeconomy.debt_kp[t] = sum(kp_id -> model[kp_id].balance.debt, all_kp)

    model.macroeconomy.debt_tot[t] = model.macroeconomy.debt_cp[t] + model.macroeconomy.debt_kp[t]

    model.macroeconomy.debt_cp_allowed[t] = Λ * sum(cp_id -> model[cp_id].curracc.S, all_cp)

    model.macroeconomy.debt_kp_allowed[t] = Λ * sum(kp_id -> model[kp_id].curracc.S, all_kp)

    if length(bankrupt_cp) > 0
        model.macroeconomy.debt_unpaid_cp[t] = sum(cp_id -> model[cp_id].balance.debt, bankrupt_cp)
    end

    if length(bankrupt_kp) > 0
        model.macroeconomy.debt_unpaid_kp[t] = sum(kp_id -> model[kp_id].balance.debt, bankrupt_kp)
    end
end


function compute_price_data!(
    # all_cp::Vector{Int},
    # all_kp::Vector{Int},
    # t::Int,
    # macroeconomy::MacroEconomy,
    model::ABM
    )

    # Compute average price, weighted by market share
    avg_p_t = sum(cp_id -> model[cp_id].p[end] * model[cp_id].f[end], model.all_cp)
    model.macroeconomy.p̄[model.t] = avg_p_t
    if model.t == 1
        model.macroeconomy.CPI[model.t] = 100
    else
        model.macroeconomy.CPI[model.t] = 100 / model.macroeconomy.p̄[begin] * avg_p_t
    end

    # Compute average price of capital goods
    avg_p_kp_t = sum(kp_id -> model[kp_id].p[end] * model[kp_id].f[end], model.all_kp)
    model.macroeconomy.p̄_kp[model.t] = avg_p_kp_t
    if model.t == 1
        model.macroeconomy.CPI_kp[model.t] = 100
    else
        model.macroeconomy.CPI_kp[model.t] = 100 / model.macroeconomy.p̄_kp[begin] * avg_p_kp_t
    end

    # Update markup rates
    model.macroeconomy.μ_cp[model.t] = mean(cp_id -> model[cp_id].μ[end], model.all_cp)
    model.macroeconomy.μ_kp[model.t] = mean(kp_id -> model[kp_id].μ[end], model.all_kp)
end


"""
Computes fraction of household that was not satisfied
"""
function compute_unsatisfied_demand(
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    all_hh::Vector{Int},
    # model.macroeconomy::MacroEconomy,
    labormarket,
    freq_per_machine::Int64,
    t::Int,
    model::ABM
    )

    model.macroeconomy.unsat_demand[t] = sum(cp_id -> model[cp_id].Dᵁ[end], all_cp) / sum(cp_id -> model[cp_id].D[end] + model[cp_id].Dᵁ[end], all_cp)

    model.macroeconomy.unspend_C[t] = 1 - sum(hh_id -> model[hh_id].C_actual, all_hh) / sum(hh_id -> model[hh_id].C, all_hh)

    model.macroeconomy.unsat_invest[t] =  1 - (sum(kp_id -> model[kp_id].Q[end], all_kp) / 
                                    (sum(cp_id -> freq_per_machine * (model[cp_id].n_mach_ordered_EI + model[cp_id].n_mach_ordered_RS), all_cp)))
    # model.macroeconomy.unsat_invest[t] = mean(cp_id -> model[cp_id].n_mach_desired_total > 0 ? (1 - (model[cp_id].n_mach_ordered_RS + model[cp_id].n_mach_ordered_EI) 
    #                                               / model[cp_id].n_mach_desired_total) : 0, all_cp)
    # total_mach_desired = 0
    # total_mach_ordered = 0
    # for cp_id in all_cp
    #     if model[cp_id].n_mach_desired_total > 0
    #         total_mach_desired += model[cp_id].n_mach_desired_total
    #         total_mach_ordered += model[cp_id].n_mach_ordered_RS + model[cp_id].n_mach_ordered_EI
    #     end
    # end
    # println(total_mach_desired, " ", total_mach_ordered)
    # model.macroeconomy.unsat_invest[t] = 1 - total_mach_ordered / total_mach_desired

    model.macroeconomy.unsat_L_demand[t] = 1 - labormarket.L_hired / labormarket.L_demanded
end


"""
Computes the GINI coefficient for wealth and income
"""
function compute_GINI(
    all_hh::Vector{Int},
    # model.macroeconomy::MacroEconomy,
    t::Int,
    model::ABM
    )

    # Compute GINI for income
    all_I = map(hh_id -> model[hh_id].total_I[end], all_hh)
    all_I_tmp = zeros(Float64, length(all_I))
    all_I_absdiff = zeros(Float64, length(all_I))

    for (i, I1) in enumerate(all_I)
        all_I_tmp .= all_I
        all_I_tmp .-= I1
        all_I_tmp .= abs.(all_I_tmp)
        all_I_absdiff[i] = sum(all_I_tmp)
    end

    model.macroeconomy.GINI_I[t] = sum(all_I_absdiff) / (2 * length(all_hh)^2 * model.macroeconomy.Ī_avg[t])

    # Compute GINI for wealth
    all_W = map(hh_id -> model[hh_id].W, all_hh)
    all_W_tmp = zeros(Float64, length(all_W))
    all_W_absdiff = zeros(Float64, length(all_W))

    for (i, W1) in enumerate(all_W)
        all_W_tmp .= all_W
        all_W_tmp .-= W1
        all_W_tmp .= abs.(all_W_tmp)
        all_W_absdiff[i] = sum(all_W_tmp)
    end

    model.macroeconomy.GINI_W[t] = sum(all_W_absdiff) / (2 * length(all_hh)^2 * model.macroeconomy.M_hh[t] / length(all_hh))
end


function compute_I_W_thresholds(
    all_hh::Vector{Int64},
    # model.macroeconomy::MacroEconomy,
    t::Int64,
    model::ABM
    )

    # Establish boundaries of middle 60%
    start_60 = round(Int64, 0.2 * length(all_hh))
    end_60 = round(Int64, 0.8 * length(all_hh))

    # Sort incomes and select income at 20th and 80th percent
    I_sorted = sort(map(hh_id -> model[hh_id].total_I, all_hh))
    model.macroeconomy.I_min[t] = I_sorted[begin]
    model.macroeconomy.I_20[t] = I_sorted[start_60]
    model.macroeconomy.I_80[t] = I_sorted[end_60]
    model.macroeconomy.I_max[t] = I_sorted[end]

    # Sort wealths and select wealth at 20th and 80th percent
    W_sorted = sort(map(hh_id -> model[hh_id].W, all_hh))
    model.macroeconomy.W_min[t] = W_sorted[begin]
    model.macroeconomy.W_20[t] = W_sorted[start_60]
    model.macroeconomy.W_80[t] = W_sorted[end_60]
    model.macroeconomy.W_max[t] = W_sorted[end]
end


function compute_savings_macro!(
    all_hh::Vector{Int64}, 
    # model.macroeconomy::MacroEconomy, 
    t::Int,
    model::ABM
    )

    all_s_emp = []
    all_s_unemp = []

    for hh_id in all_hh
        if model[hh_id].employed
            push!(all_s_emp, model[hh_id].s)
        else
            push!(all_s_unemp, model[hh_id].s)
        end
    end

    model.macroeconomy.s̄_emp[t] = length(all_s_emp) > 1 ? mean(all_s_emp) : NaN
    model.macroeconomy.s̄_unemp[t] = length(all_s_unemp) > 1 ? mean(all_s_unemp) : NaN
end


function compute_α_W_quantiles(
    all_hh::Vector{Int64},
    # model.macroeconomy::MacroEconomy, 
    t::Int64,
    model::ABM
)

    all_α = map(hh_id -> model[hh_id].α, all_hh)
    model.macroeconomy.α_W_quantiles[:, t] .= Statistics.quantile(all_α, [0.1, 0.25, 0.5, 0.75, 0.9])
end