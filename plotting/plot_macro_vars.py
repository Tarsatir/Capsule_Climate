import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars(df):
    """
    Plots macro statistics
    """
    
    fig, ax = plt.subplots(8,2, figsize=(10,19), sharex=True)

    # ax[0,0].plot(range(len(df.GDP)), df.GDP, label='total GDP')
    # ax[0,0].plot(range(len(df.GDP_I)), df.GDP_I, label='income share')
    # ax[0,0].plot(range(len(df.GDP_cp)), df.GDP_cp, label='cp share')
    # ax[0,0].plot(range(len(df.GDP_kp)), df.GDP_kp, label='kp share')
    # ax[0,0].plot(range(len(df.GDP)), df.Exp_UB, label='UB exp')
    # ax[0,0].set_title("GDP")
    # ax[0,0].legend()

    # Plot real GDP
    ax[0,0].plot(range(len(df.GDP)), 100 * df.GDP / df.CPI, label='total GDP')
    ax[0,0].plot(range(len(df.GDP_I)), 100 * df.GDP_I / df.CPI, label='income share')
    ax[0,0].plot(range(len(df.GDP_cp)), 100 * df.GDP_cp / df.CPI, label='cp share')
    ax[0,0].plot(range(len(df.GDP_kp)), 100 * df.GDP_kp / df.CPI, label='kp share')
    ax[0,0].plot(range(len(df.GDP)), 100 * df.Exp_UB / df.CPI, label='UB exp')
    ax[0,0].set_title("GDP")
    ax[0,0].legend()

    # Plot unemployment rate
    ax[0,1].plot(range(len(df.UR)), df.UR)
    ax[0,1].set_title("Unemployment rate")

    # Plot savings rate of households
    ax[1,0].plot(range(len(df.s_avg)), df.s_avg, color='red')
    ax[1,0].fill_between(range(len(df.s_avg)), df.s_avg + df.s_std, df.s_avg - df.s_std, color='red', alpha=0.4)
    ax[1,0].set_title("Savings rate")

    # Plot wage levels
    real_w_avg = 100 * df.w_avg / df.CPI
    ax[1,1].plot(range(len(real_w_avg)), real_w_avg, color='green', label='$\\bar{w}$')
    # ax[1,1].fill_between(range(len(df.w_avg)), df.w_avg + df.w_std, df.w_avg - df.w_std, 
    #                  color='green', alpha=0.4)

    real_wr_avg = 100 * df.wr_avg / df.CPI
    ax[1,1].plot(range(len(real_wr_avg)), real_wr_avg, color='red', label='$w^r$')
    # ax[1,1].fill_between(range(len(df.wr_avg)), df.wr_avg + df.wr_std, df.wr_avg - df.wr_std, 
    #                  color='red', alpha=0.4)

    real_ws_avg = 100 * df.ws_avg / df.CPI
    ax[1,1].plot(range(len(real_ws_avg)), real_ws_avg, color='blue', label='$w^s$')
    # ax[1,1].fill_between(range(len(df.ws_avg)), df.ws_avg + df.ws_std, df.ws_avg - df.ws_std, 
    #                  color='blue', alpha=0.4)
    # ax[1,1].plot(range(len(df.wo_max_avg)), df.wo_max_avg, color='orange', label='$w^o_\max$')
    ax[1,1].set_title('Real wage level')
    ax[1,1].legend()

    ax[2,0].plot(range(len(df.dL_avg)), df.dL_avg, color='red', label='all')
    ax[2,0].fill_between(range(len(df.dL_avg)), df.dL_avg + df.dL_std, df.dL_avg - df.dL_std, color='red', alpha=0.4)
    ax[2,0].plot(range(len(df.dL_cp_avg)), df.dL_cp_avg, color='green', label='cp')
    ax[2,0].plot(range(len(df.dL_kp_avg)), df.dL_kp_avg, color='orange', label='kp')
    ax[2,0].set_title('$\Delta L^d$')
    ax[2,0].legend()

    real_I_avg = 100 * df.I_avg / df.CPI
    real_I_std = 100 * df.I_std / df.I_std

    ax[2,1].plot(range(len(real_I_avg)), real_I_avg, color='purple')
    ax[2,1].fill_between(range(len(real_I_avg)), 
                         real_I_avg + real_I_std, 
                         real_I_avg - real_I_std, 
                         color='purple', 
                         alpha=0.4)
    ax[2,1].set_title('Real income of households')

    ax[3,0].plot(range(len(df.M)), df.M, label='total', zorder=10, linestyle='dashed')
    ax[3,0].plot(range(len(df.M)), df.M_hh, label='hh')
    ax[3,0].plot(range(len(df.M)), df.M_cp, label='cp')
    ax[3,0].plot(range(len(df.M)), df.M_kp, label='kp')
    ax[3,0].plot(range(len(df.M)), df.M_gov, label='gov')
    ax[3,0].set_title('Money supply')
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

    ax[4,0].plot(range(len(df.EI_avg)), 100 * df.EI_avg / df.CPI_kp, label='EI')
    ax[4,0].plot(range(len(df.RS_avg)), 100 * df.RS_avg / df.CPI_kp, label='RS')
    ax[4,0].set_title('Real investments')
    ax[4,0].legend()

    ax[4,1].plot(range(len(df.avg_pi)), df.avg_pi, label='$\\bar{\pi}$')
    ax[4,1].plot(range(len(df.avg_A)), df.avg_A, label='$\\bar{A}$')
    ax[4,1].plot(range(len(df.avg_B)), df.avg_B, label='$\\bar{B}$')
    ax[4,1].set_title('Productivity')
    ax[4,1].legend()

    ax[5,0].plot(range(len(df.avg_Q_bp)), df.avg_Q_bp, label='bp', color='blue')
    ax[5,0].plot(range(len(df.avg_Q_bp)), df.avg_n_machines_bp, 
                 label='bp n machines', color='blue', linestyle='dashed')
    ax[5,0].plot(range(len(df.avg_Q_lp)), df.avg_Q_lp, label='lp', color='red')
    ax[5,0].plot(range(len(df.avg_Q_lp)), df.avg_n_machines_lp, 
                 label='lp n machines', color='red', linestyle='dashed')
    ax[5,0].plot(range(len(df.avg_Q_kp)), df.avg_Q_kp, label='kp')
    ax[5,0].set_title('Average production quantity')
    ax[5,0].legend()

    ax[5,1].plot(range(len(df.CPI)), df.CPI, label='cp')
    ax[5,1].plot(range(len(df.CPI_kp)), df.CPI_kp, label='kp')
    ax[5,1].set_title('CPI')
    ax[5,1].legend()

    ax[6,0].plot((range(len(df.bankrupt_bp))), df.bankrupt_bp, label='bp')
    ax[6,0].plot((range(len(df.bankrupt_lp))), df.bankrupt_lp, label='lp')
    ax[6,0].plot((range(len(df.bankrupt_kp))), df.bankrupt_kp, label='kp')
    ax[6,0].legend()
    ax[6,0].set_title('Bankrupty rate')

    ax[6,1].plot(range(len(df.mu_bp)), df.mu_bp, label='bp')
    ax[6,1].plot(range(len(df.mu_lp)), df.mu_lp, label='lp')
    ax[6,1].plot(range(len(df.mu_kp)), df.mu_kp, label='kp')
    ax[6,1].legend()
    ax[6,1].set_title('Markup rates $\mu$')

    ax[7,0].plot(range(len(df.unsat_demand)), df.unsat_demand)
    ax[7,0].set_title('Unsatisfied demand')

    ax[7,1].plot(range(len(df.cu)), df.cu)
    ax[7,1].set_title('Capital utilization ratio')

    plt.tight_layout()
    plt.savefig('plots/macro_ts.png', bbox_inches='tight')


def plot_cons_vars(df):
    """
    Plots consumption figures
    """

    fig, ax = plt.subplots(2, 1, figsize=(6,4))

    # Plot real GDP growth rates
    real_GDP = 100 * df.GDP.to_numpy() / df.CPI.to_numpy()
    delta_GDP = 100 * (real_GDP[1:] - real_GDP[:-1]) / real_GDP[:-1]

    ax[0].hlines(0, 0, max(range(len(delta_GDP))), 
                 linestyle='dashed', 
                 color='black')
    ax[0].set_title('Changes in real GDP')
    ax[0].fill_between(range(len(delta_GDP)), 
                       [max(i, 0) for i in delta_GDP], 
                       [0 for _ in delta_GDP], 
                       color='green')
    ax[0].fill_between(range(len(delta_GDP)), 
                       [min(i, 0) for i in delta_GDP], 
                       [0 for _ in delta_GDP], 
                       color='red')

    ax[0].set_ylim(-7.5,7.5)


    # Plot consumption growth rates
    C_t = 100 * df.C.to_numpy() / df.CPI.to_numpy()
    delta_C = 100 * (C_t[1:] - C_t[:-1]) / C_t[:-1]

    ax[1].hlines(0, 0, max(range(len(delta_C))), linestyle='dashed', color='black')
    ax[1].set_title('Changes in real consumption')
    ax[1].fill_between(range(len(delta_C)), 
                       [max(i, 0) for i in delta_C], 
                       [0 for _ in delta_C], 
                       color='green')
    ax[1].fill_between(range(len(delta_C)), 
                       [min(i, 0) for i in delta_C], 
                       [0 for _ in delta_C], 
                       color='red')
    ax[1].set_ylim(-7.5,7.5)



    plt.tight_layout()
    plt.savefig('plots/consumption.png')


def plot_income_dist():
    """
    Plots hh income, wage and wealth distributions at end of simulation
    """

    df = pd.read_csv('../results/result_data/final_income_dists.csv')

    fig, ax = plt.subplots(3, 2, figsize=(8,8))

    ax[0,0].hist(df.all_I, bins=100)
    ax[0,0].set_title('Income ($I_{i,t}$)')

    ax[0,1].hist(df.all_I, bins=300)
    ax[0,1].set_title('Income ($I_{i,t}$), log-log')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')

    ax[1,0].hist(df.all_w, bins=100)
    ax[1,0].set_title('Wages ($w_{i,t}$)')

    ax[1,1].hist(df.all_w, bins=100)
    ax[1,1].set_title('Wages ($w_{i,t}$), log-log')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')

    ax[2,0].hist(df.all_W, bins=100)
    ax[2,0].set_title('Wealth ($W_{i,t}$)')

    ax[2,1].hist(df.all_W, bins=100)
    ax[2,1].set_title('Wealth ($W_{i,t}$), log-log')
    ax[2,1].set_xscale('log')
    ax[2,1].set_yscale('log')

    plt.tight_layout()
    plt.savefig('plots/final_income_dist.png')


def plot_sales_dist():
    """
    Plots cp and kp sales and profit distributions at end of simulation.
    """
    
    df_cp = pd.read_csv('../results/result_data/final_profit_dists_cp.csv')
    df_kp = pd.read_csv('../results/result_data/final_profit_dists_kp.csv')

    fig, ax = plt.subplots(5, 2, figsize=(8,12))

    ax[0,0].hist(df_cp.all_S_bp, bins=30)
    ax[0,0].set_title('$S$ of bp')
    
    ax[0,1].hist(df_cp.all_profit_bp, bins=30)
    ax[0,1].set_title('$\Pi$ of bp')

    ax[1,0].hist(df_cp.all_S_lp, bins=30)
    ax[1,0].set_title('$S$ of lp')
    
    ax[1,1].hist(df_cp.all_profit_lp, bins=30)
    ax[1,1].set_title('$\Pi$ of lp')

    ax[2,0].hist(df_kp.all_S_kp, bins=30)
    ax[2,0].set_title('$S$ of kp')
    
    ax[2,1].hist(df_kp.all_profit_kp, bins=30)
    ax[2,1].set_title('$\Pi$ of kp')

    ax[3,0].hist(df_cp.all_f_bp, bins=30)
    ax[3,0].set_title('$f$ of bp')
    ax[3,0].set_xlim(0, max(df_cp.all_f_bp))

    ax[3,1].hist(df_cp.all_f_lp, bins=30)
    ax[3,1].set_title('$f$ of lp')
    ax[3,1].set_xlim(0, max(df_cp.all_f_lp))

    ax[4,0].hist(df_kp.all_f_kp, bins=30)
    ax[4,0].set_title('$f$ of kp')
    ax[4,0].set_xlim(0, max(df_kp.all_f_kp))

    plt.tight_layout()
    plt.savefig('plots/final_dist_profit.png')

    
if __name__=="__main__":

    df_macro = pd.read_csv('../results/result_data/first.csv')

    plot_macro_vars(df_macro)
    plot_cons_vars(df_macro)

    plot_income_dist()
    plot_sales_dist()
