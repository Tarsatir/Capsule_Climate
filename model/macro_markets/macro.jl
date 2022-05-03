
@Base.kwdef mutable struct MacroEconomy
    T::Int=T                                                # number of timesteps
    GDP::Vector{Float64} = zeros(Float64, T)                # GDP over time
    GDP_I::Vector{Float64} = zeros(Float64, T)              # income share of GDP over time
    GDP_Π_cp::Vector{Float64} = zeros(Float64, T)           # profit share of GDP of cp over time
    GDP_Π_kp::Vector{Float64} = zeros(Float64, T)           # profit share of GDP of kp over time

    GDP_growth :: Vector{Float64} = zeros(Float64, T)       # GDP growth rates over time

    p̄::Vector{Float64} = zeros(Float64, T)                  # average price of cp goods over time
    p̄_kp::Vector{Float64} = zeros(Float64, T)               # average price of kp goods over time
    μ_cp::Vector{Float64} = zeros(Float64, T)               # average μ for cp
    μ_kp::Vector{Float64} = zeros(Float64, T)               # average μ for kp
    CPI :: Vector{Float64} = zeros(Float64, T)              # price index consumer goods over time
    CPI_kp :: Vector{Float64} = zeros(Float64, T)           # price index capital goods over time

    C::Vector{Float64} = zeros(Float64, T)                  # aggregate consumption over time
    unsat_demand::Vector{Float64} = zeros(Float64, T)       # average ratio of unsatisfied demand
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
    Ī_std::Vector{Float64} = zeros(Float64, T)              # std of income over time
    B̄_avg::Vector{Float64} = zeros(Float64, T)              # average of budget over time
    B̄_std::Vector{Float64} = zeros(Float64, T)              # std of budget over time
    U::Vector{Float64} = zeros(Float64, T)                  # unemployment over time
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
    EI_avg::Vector{Float64}  = zeros(Float64, T)            # average expansion investment
    n_mach_EI_avg::Vector{Float64} = zeros(Float64, T)      # average amount of ordered machines for EI
    RS_avg::Vector{Float64} = zeros(Float64, T)             # average replacement investment
    n_mach_RS_avg::Vector{Float64} = zeros(Float64, T)      # average amounf of ordered machines for RS

    # Productivity
    avg_π_LP::Vector{Float64} = zeros(Float64, T)           # average labor productivity cp
    avg_π_EE::Vector{Float64} = zeros(Float64, T)           # average productivity per energy unit cp

    avg_A_LP::Vector{Float64} = zeros(Float64, T)           # average A_LP at kp
    avg_A_EE::Vector{Float64} = zeros(Float64, T)           # average A_EE at kp
    avg_A_EF::Vector{Float64} = zeros(Float64, T)           # average A_EF at kp

    avg_B_LP::Vector{Float64} = zeros(Float64, T)           # average B_LP at kp
    avg_B_EE::Vector{Float64} = zeros(Float64, T)           # average B_EE at kp
    avg_B_EF::Vector{Float64} = zeros(Float64, T)           # average B_EF at kp

    # Production
    avg_Q_cp::Vector{Float64} = zeros(Float64, T)           # average production of cp
    avg_Q_kp::Vector{Float64} = zeros(Float64, T)           # average production of kp

    # Bankrupties
    bankrupt_cp::Vector{Float64} = zeros(Float64, T)        # fraction of cp that went bankrupt
    bankrupt_kp::Vector{Float64} = zeros(Float64, T)        # fraction of kp that went bankrupt

    cu::Vector{Float64} = zeros(Float64, T)                 # average rate of capital utilization
    avg_n_machines_cp::Vector{Float64} = zeros(Float64, T)  # average number of machines cp

    GINI_I::Vector{Float64} = zeros(Float64, T)             # Gini coefficient for income
    GINI_W::Vector{Float64} = zeros(Float64, T)             # Gini coefficient for wealth
end


"""
Updates macro stats after each time step
"""
function update_macro_timeseries(
    macro_struct::MacroEconomy,
    t::Int, 
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    all_kp::Vector{Int},
    all_p::Vector{Int},
    ep,
    bankrupt_cp::Vector{Int},
    bankrupt_kp::Vector{Int},
    labormarket_struct, 
    gov_struct::Government,
    indexfund_struct::IndexFund,
    global_param::GlobalParam,
    model::ABM
    )

    # Update GDP
    compute_GDP!(all_hh, all_cp, all_kp, macro_struct, t, model)

    # Update CPI
    compute_price_data!(all_cp, all_kp, t, macro_struct, model)

    # Update average labor demand
    macro_struct.ΔL̄_avg[t] = mean(p_id -> model[p_id].ΔLᵈ, all_p)
    # macro_struct.ΔL̄_std[t] = std(p_id -> model[p_id].ΔLᵈ, all_p)
    macro_struct.ΔL̄_cp_avg[t] = mean(cp_id -> model[cp_id].ΔLᵈ, all_cp)
    # macro_struct.ΔL̄_cp_std[t] = std(cp_id -> model[cp_id].ΔLᵈ, all_cp)
    macro_struct.ΔL̄_kp_avg[t] = mean(kp_id -> model[kp_id].ΔLᵈ, all_kp)

    # Consumption 
    macro_struct.C[t] = sum(hh_id->model[hh_id].C[end], all_hh)

    # Update unemployment rate and unemployment benefits expenditure
    macro_struct.U[t] = labormarket_struct.E
    macro_struct.switch_rate[t] = labormarket_struct.switch_rate
    macro_struct.Exp_UB[t] = gov_struct.curracc.Exp_UB[t]

    # Compute total amount in system
    compute_M!(all_hh, all_cp, all_kp, ep, gov_struct, indexfund_struct, 
               macro_struct, t, model)

    # Compute average savings rates
    compute_savings_macro!(all_hh, macro_struct, t, model)

    # Wage and income statistics
    update_wage_stats!(all_hh, all_p, macro_struct, t, model)

    macro_struct.Ī_avg[t] = mean(hh_id -> model[hh_id].I[end], all_hh)
    # macro_struct.Ī_std[t] = std(hh_id -> model[hh_id].I[end], all_hh)

    update_debt!(all_cp, all_kp, bankrupt_cp, bankrupt_kp, global_param.Λ, macro_struct, t, model)

    # Investment
    macro_struct.EI_avg[t] = mean(cp_id -> model[cp_id].EIᵈ, all_cp)
    macro_struct.n_mach_EI_avg[t] = mean(cp_id -> model[cp_id].n_mach_ordered_EI, all_cp)
    macro_struct.RS_avg[t] = mean(cp_id -> model[cp_id].RSᵈ, all_cp)
    macro_struct.n_mach_RS_avg[t] = mean(cp_id -> model[cp_id].n_mach_ordered_RS, all_cp)

    # Productivity
    macro_struct.avg_π_LP[t] = mean(cp_id -> model[cp_id].π_LP, all_cp)
    macro_struct.avg_π_EE[t] = mean(cp_id -> model[cp_id].π_EE, all_cp)

    macro_struct.avg_A_LP[t] = mean(kp_id -> model[kp_id].A_LP[end], all_kp)
    macro_struct.avg_A_EE[t] = mean(kp_id -> model[kp_id].A_EE[end], all_kp)
    macro_struct.avg_A_EF[t] = mean(kp_id -> model[kp_id].A_EF[end], all_kp)

    macro_struct.avg_B_LP[t] = mean(kp_id -> model[kp_id].B_LP[end], all_kp)
    macro_struct.avg_B_EE[t] = mean(kp_id -> model[kp_id].B_EE[end], all_kp)
    macro_struct.avg_B_EF[t] = mean(kp_id -> model[kp_id].B_EF[end], all_kp)

    # Production quantity
    macro_struct.avg_Q_cp[t] = mean(cp_id -> model[cp_id].Q[end], all_cp)
    macro_struct.avg_Q_kp[t] = mean(kp_id -> model[kp_id].Q[end], all_kp)

    compute_bankrupties(all_cp, all_kp, bankrupt_cp, bankrupt_kp, macro_struct, t)

    compute_unsatisfied_demand(all_cp, all_hh, macro_struct, t, model)

    macro_struct.avg_N_goods[t] = mean(cp_id -> model[cp_id].N_goods, all_cp)

    # Mean rate of capital utilization
    macro_struct.cu[t] = mean(cp_id -> model[cp_id].n_machines > 0 ? model[cp_id].cu : 0.5, all_cp)

    # Average number of machines
    macro_struct.avg_n_machines_cp[t] = mean(cp_id -> model[cp_id].n_machines, all_cp)

    # Compute GINI coefficients
    compute_GINI(all_hh, macro_struct, t, model)

end


"""
Computes GDP based on income of separate sectors, computes partial incomes of sectors
"""
function compute_GDP!(
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    # Household income
    total_I = sum(hh_id -> model[hh_id].Iᵀ[end], all_hh)
    macro_struct.GDP_I[t] = total_I

    # cp profits
    total_Π_cp = sum(cp_id -> model[cp_id].Π[end], all_cp)
    macro_struct.GDP_Π_cp[t] = total_Π_cp
    
    # kp profits
    total_Π_kp = sum(kp_id -> model[kp_id].Π[end], all_kp)
    macro_struct.GDP_Π_kp[t] = total_Π_kp

    # total GDP
    macro_struct.GDP[t] = total_I + total_Π_cp + total_Π_kp

    # GDP growth
    if t > 1
        macro_struct.GDP_growth[t] = (macro_struct.GDP[t] - macro_struct.GDP[t-1]) / macro_struct.GDP[t-1]
    end
end


"""
Computes the ratios of bankrupt bp, lp and kp.
"""
function compute_bankrupties(
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    bankrupt_cp::Vector{Int},
    bankrupt_kp::Vector{Int},
    macro_struct::MacroEconomy,
    t::Int
    )

    macro_struct.bankrupt_cp[t] = length(bankrupt_cp) / length(all_cp)

    macro_struct.bankrupt_kp[t] = length(bankrupt_kp) / length(all_kp)

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
    gov_struct::Government,
    indexfund_struct::IndexFund,
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    # Wealth of households
    macro_struct.M_hh[t] = sum(hh_id -> model[hh_id].W, all_hh)

    # Liquid assets of cp_id
    macro_struct.M_cp[t] = sum(cp_id -> model[cp_id].balance.NW, all_cp)

    # Liquid assets of kp
    macro_struct.M_kp[t] = sum(kp_id -> model[kp_id].balance.NW, all_kp)

    # Liquid assets of ep
    macro_struct.M_ep[t] = ep.NWₑ[t]

    # Money owned by government
    macro_struct.M_gov[t] = gov_struct.MS

    # Money in investment fund
    macro_struct.M_if[t] = indexfund_struct.Assets

    # Total amount of money stocks
    macro_struct.M[t] = (macro_struct.M_hh[t] + macro_struct.M_cp[t] + macro_struct.M_kp[t] +
                         macro_struct.M_ep[t] + macro_struct.M_gov[t] + macro_struct.M_if[t])
end


"""
Computes average wage statistics
"""
function update_wage_stats!(
    all_hh::Vector{Int},
    all_p::Vector{Int},
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    macro_struct.w̄_avg[t] = mean(p_id -> model[p_id].w̄[end], all_p)

    macro_struct.wʳ_avg[t] = mean(hh_id -> model[hh_id].wʳ, all_hh)

    macro_struct.wˢ_avg[t] = mean(hh_id -> model[hh_id].wˢ, all_hh)

    macro_struct.wᴼ_max_mean[t] = mean(p_id -> model[p_id].wᴼ_max, all_p)
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
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    macro_struct.debt_cp[t] = sum(cp_id -> model[cp_id].balance.debt, all_cp)

    macro_struct.debt_kp[t] = sum(kp_id -> model[kp_id].balance.debt, all_kp)

    macro_struct.debt_tot[t] = macro_struct.debt_cp[t] + macro_struct.debt_kp[t]

    macro_struct.debt_cp_allowed[t] = Λ * sum(cp_id -> model[cp_id].curracc.S, all_cp)

    macro_struct.debt_kp_allowed[t] = Λ * sum(kp_id -> model[kp_id].curracc.S, all_kp)

    if length(bankrupt_cp) > 0
        macro_struct.debt_unpaid_cp[t] = sum(cp_id -> model[cp_id].balance.debt, bankrupt_cp)
    end

    if length(bankrupt_kp) > 0
        macro_struct.debt_unpaid_kp[t] = sum(kp_id -> model[kp_id].balance.debt, bankrupt_kp)
    end
end


function compute_price_data!(
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    t::Int,
    macro_struct::MacroEconomy,
    model::ABM
    )

    # Compute average price, weighted by market share
    avg_p_t = sum(cp_id -> model[cp_id].p[end] * model[cp_id].f[end], all_cp)
    macro_struct.p̄[t] = avg_p_t

    if t == 1
        macro_struct.CPI[t] = 100
    else
        macro_struct.CPI[t] = 100 / macro_struct.p̄[1] * avg_p_t
    end

    # Compute average price of capital goods
    avg_p_kp_t = sum(kp_id -> model[kp_id].p[end] * model[kp_id].f[end], all_kp)
    macro_struct.p̄_kp[t] = avg_p_kp_t
    if t == 1
        macro_struct.CPI_kp[t] = 100
    else
        macro_struct.CPI_kp[t] = 100 / macro_struct.p̄_kp[1] * avg_p_kp_t
    end

    # Update markup rates
    macro_struct.μ_cp[t] = mean(cp_id -> model[cp_id].μ[end], all_cp)
    macro_struct.μ_kp[t] = mean(kp_id -> model[kp_id].μ[end], all_kp)
end


"""
Computes fraction of household that was not satisfied
"""
function compute_unsatisfied_demand(
    all_cp::Vector{Int},
    all_hh::Vector{Int},
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    # mean_unsat_dem = 0.0

    # for hh_id in all_hh
    #     if length(model[hh_id].unsat_dem) > 0 
    #         mean_unsat_dem += mean(ud->ud[2], model[hh_id].unsat_dem)
    #     end
    # end

    # macro_struct.unsat_demand[t] = mean_unsat_dem / length(all_hh)
    macro_struct.unsat_demand[t] = sum(cp_id -> model[cp_id].Dᵁ, all_cp) / sum(cp_id -> model[cp_id].D[end] + model[cp_id].Dᵁ, all_cp)
end


"""
Computes the GINI coefficient for wealth and income
"""
function compute_GINI(
    all_hh::Vector{Int},
    macro_struct::MacroEconomy,
    t::Int,
    model::ABM
    )

    # Compute GINI for income
    all_I = map(hh_id -> model[hh_id].Iᵀ[end], all_hh)
    all_I_absdiff = zeros(length(all_I))

    for (i, I1) in enumerate(all_I)
        for I2 in all_I
            @inbounds all_I_absdiff[i] += abs(I1 - I2)
        end
    end
    macro_struct.GINI_I[t] = sum(all_I_absdiff) / (2 * length(all_hh)^2 * macro_struct.Ī_avg[t])

    # Compute GINI for wealth
    all_W = map(hh_id -> model[hh_id].W, all_hh)
    all_W_absdiff = zeros(length(all_W))

    for (i, W1) in enumerate(all_W)
        for W2 in all_W
            @inbounds all_W_absdiff[i] += abs(W1 - W2)
        end
    end
    macro_struct.GINI_W[t] = sum(all_W_absdiff) / (2 * length(all_hh)^2 * macro_struct.M_hh[t] / length(all_hh))
end


function compute_savings_macro!(
    all_hh::Vector{Int}, 
    macro_struct::MacroEconomy, 
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

    macro_struct.s̄_emp[t] = length(all_s_emp) > 1 ? mean(all_s_emp) : NaN
    macro_struct.s̄_unemp[t] = length(all_s_unemp) > 1 ? mean(all_s_unemp) : NaN
end