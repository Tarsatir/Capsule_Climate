from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars(df):
    """
    Plots macro statistics
    """
    
    fig, ax = plt.subplots(8,2, figsize=(10,25), sharex=True)

    T = range(len(df.GDP))

    # ax[0,0].plot(range(len(df.GDP)), df.GDP, label='total GDP')
    # ax[0,0].plot(range(len(df.GDP_I)), df.GDP_I, label='income share')
    # ax[0,0].plot(range(len(df.GDP_cp)), df.GDP_cp, label='cp share')
    # ax[0,0].plot(range(len(df.GDP_kp)), df.GDP_kp, label='kp share')
    # ax[0,0].plot(range(len(df.GDP)), df.Exp_UB, label='UB exp')
    # ax[0,0].set_title("GDP")
    # ax[0,0].legend()

    # Plot real GDP
    ax[0,0].plot(T, 100 * df.GDP / df.CPI, label='total GDP')
    ax[0,0].plot(T, 100 * df.GDP_I / df.CPI, label='income share')
    ax[0,0].plot(T, 100 * df.GDP_cp / df.CPI, label='cp share')
    ax[0,0].plot(T, 100 * df.GDP_kp / df.CPI, label='kp share')
    ax[0,0].plot(T, 100 * df.Exp_UB / df.CPI, label='UB exp')
    ax[0,0].set_title("GDP")
    ax[0,0].legend()

    # ax[0,0].plot(T, df.GDP, label='total GDP')
    # ax[0,0].plot(T, df.GDP_I, label='income share')
    # ax[0,0].plot(T, df.GDP_cp, label='cp share')
    # ax[0,0].plot(T, df.GDP_kp, label='kp share')
    # ax[0,0].plot(T, df.Exp_UB, label='UB exp')
    # ax[0,0].set_title("GDP")
    # ax[0,0].legend()

    # Plot unemployment rate
    ax[0,1].plot(range(len(df.UR)), df.UR, label='unemployment rate')
    ax[0,1].plot(range(len(df.switch_rate)), df.switch_rate, label='switching rate')
    ax[0,1].set_title("Unemployment rate")
    ax[0,1].set_ylim(0,1)
    ax[0,1].legend()

    ax[1,0].plot(range(len(df.M)), df.M - df.debt_tot, 
                 label='total', zorder=5, linestyle='dashed')
    ax[1,0].plot(range(len(df.M)), df.M_hh, label='hh')
    ax[1,0].plot(range(len(df.M)), df.M_cp, label='cp')
    ax[1,0].plot(range(len(df.M)), df.M_kp, label='kp')
    ax[1,0].plot(range(len(df.M)), df.M_ep, label='ep')
    ax[1,0].plot(range(len(df.M)), df.M_gov, label='gov')
    ax[1,0].plot(range(len(df.M)), df.debt_tot, label='total debts')
    ax[1,0].hlines(df.M[0], 0, len(df.M), linestyle='dotted', alpha=0.5, color='black')
    # ax[3,0].plot(range(len(df.M)), -np.cumsum(df.debt_unpaid_cp.to_numpy()), label='unpaid cp debt')
    # ax[3,0].plot(range(len(df.M)), -np.cumsum(df.debt_unpaid_kp.to_numpy()), label='unpaid kp debt')
    ax[1,0].plot(range(len(df.M)), df.M_if, label='if')
    ax[1,0].set_title('Money supply')
    ax[1,0].legend()

    # Plot wage levels
    real_w_avg = 100 * df.w_avg / df.CPI
    # real_w_avg = 100 * df.w_avg
    ax[1,1].plot(range(len(real_w_avg)), real_w_avg, color='green', label='real $\\bar{w}$')
    ax[1,1].plot(range(len(df.w_avg)), df.w_avg, color='blue', label='$\\bar{w}$')
    ax[1,1].set_title('Real wage level')
    ax[1,1].legend()

    ax[2,0].hlines(0, 0, T[-1], linestyle='dashed', color='black')
    ax[2,0].plot(range(len(df.dL_avg)), df.dL_avg, color='red', label='all')
    ax[2,0].fill_between(range(len(df.dL_avg)), df.dL_avg + df.dL_std, df.dL_avg - df.dL_std, color='red', alpha=0.4)
    ax[2,0].plot(range(len(df.dL_cp_avg)), df.dL_cp_avg, color='green', label='cp')
    ax[2,0].plot(range(len(df.dL_kp_avg)), df.dL_kp_avg, color='orange', label='kp')
    ax[2,0].set_title('$\Delta L^d$')
    ax[2,0].legend()

    real_I_avg = 100 * df.I_avg / df.CPI
    real_I_labor = 100 * df.I_labor_avg / df.CPI
    real_I_capital = 100 * df.I_capital_avg / df.CPI
    real_I_UB = 100 * df.I_UB_avg / df.CPI
    real_I_socben = 100 * df.I_socben_avg / df.CPI

    ax[2,1].plot(T, real_I_avg, color='purple', label='total')
    ax[2,1].plot(T, real_I_labor, color='blue', label='labor')
    ax[2,1].plot(T, real_I_capital, color='red', label='capital')
    ax[2,1].plot(T, real_I_UB, color='green', label='UB')
    ax[2,1].plot(T, real_I_socben, color='orange', label='socben')
    ax[2,1].hlines(0, max(T), 0, linestyle='dashed', color='black')
    ax[2,1].legend()
    ax[2,1].set_title('Real income of households')
    

    # Plot savings rate of households
    ax[3,0].hlines(0, 0, T[-1], linestyle='dashed', color='black')
    ax[3,0].plot(T, df.s_emp, color='red', label='employed')
    ax[3,0].plot(T, df.s_unemp, color='blue', label='unemployed')
    ax[3,0].set_title("Savings rate")
    ax[3,0].legend()


    ax[3,1].plot(range(len(df.debt_tot)), df.debt_tot, label='total')
    ax[3,1].plot(range(len(df.debt_cp)), df.debt_cp, label='cp', color='green')
    ax[3,1].plot(range(len(df.debt_cp_allowed)), df.debt_cp_allowed, 
                 label='cp allowed', color='green', linestyle='dashed', alpha=0.5)
    ax[3,1].plot(range(len(df.debt_kp)), df.debt_kp, label='kp', color='red')
    ax[3,1].plot(range(len(df.debt_kp_allowed)), df.debt_kp_allowed, 
                label='kp allowed', color='red', linestyle='dashed', alpha=0.5)
    ax[3,1].set_title('Debt levels')
    ax[3,1].legend()

    # ax[4,0].plot(range(len(df.EI_avg)), 100 * df.EI_avg / df.CPI_kp, label='EI')
    # ax[4,0].plot(range(len(df.RS_avg)), 100 * df.RS_avg / df.CPI_kp, label='RS')
    # ax[4,0].set_title('Real investments')
    ax[4,0].plot(range(len(df.n_mach_EI)), df.n_mach_EI, label='n EI')
    ax[4,0].plot(range(len(df.n_mach_RS)), df.n_mach_RS, label='n RS')
    ax[4,0].legend()

    ax[4,1].plot(range(len(df.avg_pi_LP)), df.avg_pi_LP, label='$\\bar{\pi}_{LP}$')
    ax[4,1].plot(range(len(df.avg_pi_EE)), df.avg_pi_EE, label='$\\bar{\pi}_{EE}$')
    ax[4,1].plot(range(len(df.avg_A_LP)), df.avg_A_LP, label='$\\bar{A}_{LP}$')
    ax[4,1].plot(range(len(df.avg_A_EE)), df.avg_A_EE, label='$\\bar{A}_{EE}$')
    ax[4,1].plot(range(len(df.avg_A_EF)), df.avg_A_EF, label='$\\bar{A}_{EF}$')
    ax[4,1].plot(range(len(df.avg_B_LP)), df.avg_B_LP, label='$\\bar{B}_{LP}}$')
    ax[4,1].plot(range(len(df.avg_B_EE)), df.avg_B_EE, label='$\\bar{B}_{EE}}$')
    ax[4,1].plot(range(len(df.avg_B_EF)), df.avg_B_EF, label='$\\bar{B}_{EF}}$')
    ax[4,1].set_title('Productivity')
    ax[4,1].legend()


    ax[5,0].plot(T, df.avg_Q_kp, label='kp Q')
    ax[5,0].plot(T, df.avg_Q_cp, label='cp Q', color='green')
    ax[5,0].plot(T, df.avg_Qs_cp, label='cp $Q^s$', color='green', linestyle='dashed')
    ax[5,0].plot(T, df.avg_Qe_cp, label='cp $Q^e$', color='green', linestyle='dotted')
    ax[5,0].plot(T, df.avg_n_machines_cp, 
                 label='cp n machines', color='blue', linestyle='dashed')
    ax[5,0].plot(T, df.avg_D_cp, label='cp $D$', color='red')
    ax[5,0].plot(T, df.avg_De_cp, label='cp $D^e$', color='red', linestyle='dashed')
    ax[5,0].plot(T, df.avg_Du_cp, label='cp $D^U$', color='red', linestyle='dotted', alpha=0.7)
    # ax[5,0].plot(range(len(df.avg_Q_lp)), df.avg_Q_lp, label='lp', color='red')
    # ax[5,0].plot(range(len(df.avg_Q_lp)), df.avg_n_machines_lp, 
    #              label='lp n machines', color='red', linestyle='dashed')
    ax[5,0].plot(T, df.avg_N_goods, label='avg $N$', color='orange')
    ax[5,0].set_title('Average production quantity')
    ax[5,0].legend()

    ax[5,1].plot(range(len(df.CPI)), df.CPI, label='cp')
    ax[5,1].plot(range(len(df.CPI_kp)), df.CPI_kp, label='kp')
    ax[5,1].set_title('CPI')
    ax[5,1].legend()


    ax[6,0].plot(T, df.bankrupt_kp, label='kp')
    ax[6,0].plot(T, df.bankrupt_cp, label='cp')
    ax[6,0].legend()
    ax[6,0].set_title('Bankrupty rate')

    ax[6,1].plot(range(len(df.mu_cp)), df.mu_cp, label='cp')
    # ax[6,1].plot(range(len(df.mu_lp)), df.mu_lp, label='lp')
    ax[6,1].plot(range(len(df.mu_kp)), df.mu_kp, label='kp')
    ax[6,1].legend()
    ax[6,1].set_title('Markup rates $\mu$')

    ax[7,0].plot(T, df.unsat_demand, label='unsatisfied D')
    ax[7,0].plot(T, df.unspend_C, label='unspend C')
    ax[7,0].plot(T, df.unsat_invest, label='unsatisfied I')
    ax[7,0].plot(T, df.unsat_L_demand, label='unsatisfied L')
    ax[7,0].plot(T, df.cu, label='cu')
    ax[7,0].set_title('Unsatisfied demand and unspend C')
    ax[7,0].set_ylim(0, 1)
    ax[7,0].legend()

    ax[7,1].plot(T, df.total_C, label="$C$")
    ax[7,1].plot(T, df.total_C_actual, label="$C$ actual")
    ax[7,1].plot(T, df.total_I, label="$I$")
    ax[7,1].plot(T, df.total_w, label="$w$")
    ax[7,1].set_title('Spending')
    ax[7,1].legend()

    plt.tight_layout()
    plt.savefig('plots/macro_ts.png', bbox_inches='tight')


def plot_cons_vars(df):
    """
    Plots consumption figures
    """

    fig = plt.figure(figsize=(6, 12))

    gs = GridSpec(6, 2, figure=fig)
    
    ax0 = fig.add_subplot(gs[0,:])
    ax1 = fig.add_subplot(gs[1,:])
    ax2 = fig.add_subplot(gs[2,:])
    ax3 = fig.add_subplot(gs[3,:])
    ax4 = fig.add_subplot(gs[4,:])
    ax5 = fig.add_subplot(gs[5,0])
    ax6 = fig.add_subplot(gs[5,1])

    
    if len(df.GDP) <= 100:
        return

    # Plot real GDP growth rates
    real_GDP = 100 * df.GDP.to_numpy()[100:] / df.CPI.to_numpy()[100:]
    delta_GDP = 100 * (real_GDP[1:] - real_GDP[:-1]) / real_GDP[:-1]

    T = np.arange(100, 100+len(real_GDP)-1)
    

    ax0.hlines(0, min(T), max(T), linestyle='dashed', color='black')
    ax0.set_title('Monhtly changes in real GDP')
    ax0.fill_between(T, [max(i, 0) for i in delta_GDP], 
                        [0 for _ in delta_GDP], color='green')
    ax0.fill_between(T, [min(i, 0) for i in delta_GDP], 
                        [0 for _ in delta_GDP], color='red')
    ax0.set_xlabel('time')
    ax0.set_ylabel('growth rate (%)')
    ax0.set_ylim(-7.5,7.5)


    # Plot consumption growth rates
    C_t = 100 * df.C.to_numpy()[100:] / df.CPI.to_numpy()[100:]
    delta_C = 100 * (C_t[1:] - C_t[:-1]) / C_t[:-1]

    ax1.hlines(0, min(T), max(T), linestyle='dashed', color='black')
    ax1.set_title('Monhtly changes in real consumption')
    ax1.fill_between(T, [max(i, 0) for i in delta_C], 
                        [0 for _ in delta_C], color='green')
    ax1.fill_between(T, [min(i, 0) for i in delta_C], 
                        [0 for _ in delta_C], color='red')
    ax1.set_ylabel('growth rate (%)')
    ax1.set_xlabel('time')
    ax1.set_ylim(-7.5,7.5)

    # Plot hh GDP rowth rates
    real_GDP_hh = 100 * df.GDP_I.to_numpy()[100:] / df.CPI.to_numpy()[100:]
    delta_GDP_hh = 100 * (real_GDP_hh[1:] - real_GDP_hh[:-1]) / real_GDP_hh[:-1]

    ax2.hlines(0, min(T), max(T), linestyle='dashed', color='black')
    ax2.set_title('Monthly changes in real GDP of hh')
    ax2.fill_between(T, [max(i, 0) for i in delta_GDP_hh], 
                        [0 for _ in delta_GDP_hh], color='green')
    ax2.fill_between(T, [min(i, 0) for i in delta_GDP_hh], 
                        [0 for _ in delta_GDP_hh], color='red')
    ax2.set_xlabel('time')
    ax2.set_ylabel('growth rate (%)')
    ax2.set_ylim(-7.5,7.5)

    # Plot cp GDP growth rates
    real_GDP_cp = 100 * df.GDP_cp.to_numpy()[100:] / df.CPI.to_numpy()[100:]
    delta_GDP_cp = 100 * (real_GDP_cp[1:] - real_GDP_cp[:-1]) / real_GDP_cp[:-1]

    ax3.hlines(0, min(T), max(T), linestyle='dashed', color='black')
    ax3.set_title('Monthly changes in real GDP of cp')
    ax3.fill_between(T, [max(i, 0) for i in delta_GDP_cp], 
                        [0 for _ in delta_GDP_cp], color='green')
    ax3.fill_between(T, [min(i, 0) for i in delta_GDP_cp], 
                        [0 for _ in delta_GDP_cp], color='red')
    ax3.set_xlabel('time')
    ax3.set_ylabel('growth rate (%)')
    ax3.set_ylim(-7.5,7.5)

    # Plot kp GDP growth rates
    real_GDP_kp = 100 * df.GDP_kp.to_numpy()[100:] / df.CPI.to_numpy()[100:]
    delta_GDP_kp = 100 * (real_GDP_kp[1:] - real_GDP_kp[:-1]) / real_GDP_kp[:-1]

    ax4.hlines(0, min(T), max(T), linestyle='dashed', color='black')
    ax4.set_title('Monthly changes in real GDP of kp')
    ax4.fill_between(T, [max(i, 0) for i in delta_GDP_kp], 
                        [0 for _ in delta_GDP_kp], color='green')
    ax4.fill_between(T, [min(i, 0) for i in delta_GDP_kp], 
                        [0 for _ in delta_GDP_kp], color='red')
    ax4.set_xlabel('time')
    ax4.set_ylabel('growth rate (%)')
    ax4.set_ylim(-7.5,7.5)

    # Compute quarterly GDP growth rates and plot
    Q_delta_GDP = 100 * (real_GDP[3:] - real_GDP[:-3]) / real_GDP[:-3]

    ax5.set_title('Quarterly GDP growth')
    ax5.hist(Q_delta_GDP, bins=100, density=True)
    ax5.set_xlabel('growth rate (%)')
    ax5.set_xlim(-7.5, 7.5)

    # Compute quarterly C growth rates and plot
    Q_delta_C = 100 * (C_t[3:] - C_t[:-3]) / C_t[:-3]
    ax6.set_title('Quarterly $C$ growth')
    ax6.hist(Q_delta_C, bins=100, density=True)
    ax6.set_xlabel('growth rate (%)')
    ax6.set_xlim(-7.5, 7.5)


    plt.tight_layout()
    plt.savefig('plots/consumption.png')


def plot_income_dist():
    """
    Plots hh income, wage and wealth distributions at end of simulation
    """

    df = pd.read_csv('../results/result_data/final_income_dists.csv')

    fig, ax = plt.subplots(4, 2, figsize=(8,10))

    ax[0,0].hist(df.all_I, bins=100)
    ax[0,0].set_title('Income ($I_{i,t}$)')
    ax[0,0].set_xlim(0, max(df.all_I))

    ax[0,1].hist(df.all_I, bins=300)
    ax[0,1].set_title('Income ($I_{i,t}$), log-log')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')

    ax[1,0].hist(df.all_w, bins=100)
    ax[1,0].set_title('Wages ($w_{i,t}$)')
    ax[1,0].set_xlim(0, max(df.all_w))

    ax[1,1].hist(df.all_w, bins=100)
    ax[1,1].set_title('Wages ($w_{i,t}$), log-log')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')

    ax[2,0].hist(df.all_W, bins=100)
    ax[2,0].set_title('Wealth ($W_{i,t}$)')
    ax[2,0].set_xlim(0, max(df.all_W))

    ax[2,1].hist(df.all_W, bins=100)
    ax[2,1].set_title('Wealth ($W_{i,t}$), log-log')
    ax[2,1].set_xscale('log')
    ax[2,1].set_yscale('log')

    ax[3,0].scatter(df.skills, df.all_I, s=0.3)
    ax[3,0].set_title('Income to skills')
    ax[3,0].set_xlabel('skill')
    ax[3,0].set_ylabel('income')

    ax[3,1].scatter(df.skills, df.all_W, s=0.3)
    ax[3,1].set_title('Wealth to skills')
    ax[3,1].set_xlabel('skill')
    ax[3,1].set_ylabel('wealth')

    plt.tight_layout()
    plt.savefig('plots/final_income_dist.png')


def plot_sales_dist():
    """
    Plots cp and kp sales and profit distributions at end of simulation.
    """
    
    df_cp = pd.read_csv('../results/result_data/final_profit_dists_cp.csv')
    df_kp = pd.read_csv('../results/result_data/final_profit_dists_kp.csv')

    fig, ax = plt.subplots(5, 2, figsize=(8,12))

    ax[0,0].hist(df_cp.all_S_cp, bins=60)
    ax[0,0].set_title('$S$ of cp')
    
    ax[1,0].hist(df_cp.all_profit_cp, bins=60)
    ax[1,0].set_title('$\Pi$ of cp')

    ax[2,0].hist(df_cp.all_f_cp, bins=60)
    ax[2,0].set_title('$f$ of cp')
    ax[2,0].set_xlim(0, max(df_cp.all_f_cp))

    ax[3,0].hist(df_cp.all_L_cp, bins=30)
    ax[3,0].set_title("$L$ of cp")

    ax[0,1].hist(df_kp.all_S_kp, bins=30)
    ax[0,1].set_title('$S$ of kp')
    
    ax[1,1].hist(df_kp.all_profit_kp, bins=30)
    ax[1,1].set_title('$\Pi$ of kp')

    ax[2,1].hist(df_kp.all_f_kp, bins=30)
    ax[2,1].set_title('$f$ of kp')
    ax[2,1].set_xlim(0, max(df_kp.all_f_kp))

    ax[3,1].hist(df_kp.all_L_kp, bins=30)
    ax[3,1].set_title("$L$ of kp")

    ax[4,0].scatter(df_cp.all_w_cp, df_cp.all_L_cp)
    ax[4,0].set_title("$w$ to $L$")

    ax[4,1].scatter(df_cp.all_p_cp, df_cp.all_L_cp)
    ax[4,1].set_title("$p$ to $L$")

    plt.tight_layout()
    plt.savefig('plots/final_dist_profit.png')


def plot_inequality(df_macro):
    """
    Plot GINI coefficients for income and wealth over time
    """

    fig, ax = plt.subplots(1, 2, figsize=(8,4))

    T = range(len(df_macro.gini_I))

    ax[0].plot(T, df_macro.gini_I, label='model output')
    ax[0].plot(T, df_macro.FGT, label='FGT index')
    ax[0].hlines(0.282, 0, len(df_macro.gini_I), linestyle='dashed', color='black', 
                 label='Netherlands (2018)')
    ax[0].set_ylim(0,1)
    ax[0].set_title("Income inequality")
    ax[0].legend()

    ax[1].plot(T, df_macro.gini_W, label='model output')
    ax[1].hlines(0.789, 0, len(df_macro.gini_W), linestyle='dashed', color='black', 
                 label='Netherlands (2018)')
    ax[1].set_title("Wealth inequality")
    ax[1].set_ylim(0,1)
    ax[1].legend()

    plt.tight_layout()
    plt.savefig('plots/inequality.png')


def plot_energy(df_climate_energy, df_macro):

    fig, ax = plt.subplots(2, 2, figsize=(8,6))

    T = range(len(df_climate_energy.energy_demand))

    # Plot energy use and capacities
    ax[0,0].plot(T, df_climate_energy.energy_demand, label='$D_{e}(t)$', color='red')
    ax[0,0].plot(T, df_climate_energy.total_capacity, label='$\\bar{Q}_e$', 
               color='blue', linestyle='dashed')
    ax[0,0].plot(T, df_climate_energy.green_capacity, label='green capacity', 
               color='green')
    ax[0,0].plot(T, df_climate_energy.total_capacity - df_climate_energy.green_capacity, 
               label='dirty capacity', color='brown')
    ax[0,0].set_title('Energy demand and consumption')
    ax[0,0].set_xlabel('Time')
    ax[0,0].set_ylabel('Units of energy')
    ax[0,0].legend()

    # Plot energy intensity
    ax[0,1].plot(T, df_climate_energy.energy_demand / (df_macro.GDP / df_macro.CPI))
    ax[0,1].set_title('Energy intensity per unit of real GDP')
    ax[0,1].set_xlabel('Time')
    ax[0,1].set_ylabel('Energy intensity')
    
    # Plot innovation spending
    ax[1,0].plot(T, df_climate_energy.RD, label='total R&D spending', color='blue', linestyle='dashed')
    ax[1,0].plot(T, df_climate_energy.IN_g, label='green R&D spending', color='green')
    ax[1,0].plot(T, df_climate_energy.IN_d, label='dirty R&D spending', color='brown')
    ax[1,0].legend()

    ax[1,1].plot(T, df_climate_energy.p_e, label='energy prices')
    ax[1,1].set_title('Energy prices')
    ax[1,1].legend()

    plt.tight_layout()
    plt.savefig('plots/energy.png')


def plot_climate(df_climate_energy, df_macro):

    fig, ax = plt.subplots(2, 2, figsize=(8,6))

    T = range(len(df_climate_energy.emissions_total))

    ax[0,0].plot(T, df_climate_energy.emissions_total, label='$c^{total}_t$')
    ax[0,0].plot(T, df_climate_energy.emissions_kp, label='$c^{kp}_t$')
    ax[0,0].plot(T, df_climate_energy.emissions_cp, label='$c^{cp}_t$')
    ax[0,0].plot(T, df_climate_energy.emissions_ep, label='$c^{ep}_t$')
    ax[0,0].set_title('CO$_2$ emissions')
    ax[0,0].set_xlabel('time')
    ax[0,0].set_ylabel('total CO$_2$ emission')
    ax[0,0].legend()

    real_GDP = 100 * df_macro.GDP / df_macro.CPI
    ax[0,1].plot(T, df_climate_energy.emissions_total / real_GDP, label='total emissions')
    ax[0,1].plot(T, df_climate_energy.emissions_kp / real_GDP, label='kp emissions')
    ax[0,1].plot(T, df_climate_energy.emissions_cp / real_GDP, label='cp emissions')
    ax[0,1].plot(T, df_climate_energy.emissions_ep / real_GDP, label='ep emissions')
    ax[0,1].set_title('CO$_2$ emissions per unit of real GDP')
    ax[0,1].set_xlabel('time')
    ax[0,1].set_ylabel('CO$_2$ / GDP')
    ax[0,1].legend()


    ax[1,0].plot(T, df_climate_energy.C_a, label='CO$_2$ in atmosphere')
    ax[1,0].plot(T, df_climate_energy.C_m, label='CO$_2$ in mixed ocean layer')
    ax[1,0].plot(T, df_climate_energy.C_d, label='CO$_2$ in deep ocean layer')
    # ax[1,0].plot(T, df_climate_energy.NPP, label='NPP$_t$')
    ax[1,0].set_title('CO$_2$ concentrations')
    ax[1,0].set_xlabel('time')
    ax[1,0].set_ylabel('Total CO$_2$ concentration')
    # ax[1,0].set_yscale('log')
    ax[1,0].legend()

    ax[1,1].plot(T, df_climate_energy.dT_m, label='$\delta T_{m,t}$')
    ax[1,1].plot(T, df_climate_energy.dT_d, label='$\delta T_{d,t}$')
    ax[1,1].set_title('Temperatures')
    ax[1,1].set_xlabel('time')
    ax[1,1].set_ylabel('Temperature anomaly')
    ax[1,1].legend()

    plt.tight_layout()
    plt.savefig('plots/climate.png')

    
if __name__=="__main__":

    df_macro = pd.read_csv('../results/result_data/first.csv')

    plot_macro_vars(df_macro)
    plot_cons_vars(df_macro)

    plot_income_dist()
    plot_inequality(df_macro)
    plot_sales_dist()

    df_climate_energy = pd.read_csv('../results/result_data/climate_and_energy.csv')
    plot_energy(df_climate_energy, df_macro)
    plot_climate(df_climate_energy, df_macro)
