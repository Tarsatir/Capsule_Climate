import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars():
    """
    Plots macro statistics
    """

    df = pd.read_csv('../results/result_data/first.csv')
    
    fig, ax = plt.subplots(7,2, figsize=(10,17), sharex=True)

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
    ax[1,1].plot(range(len(df.w_avg)), df.w_avg, color='green', label='$\\bar{w}$')
    # ax[1,1].fill_between(range(len(df.w_avg)), df.w_avg + df.w_std, df.w_avg - df.w_std, 
    #                  color='green', alpha=0.4)
    ax[1,1].plot(range(len(df.wr_avg)), df.wr_avg, color='red', label='$w^r$')
    # ax[1,1].fill_between(range(len(df.wr_avg)), df.wr_avg + df.wr_std, df.wr_avg - df.wr_std, 
    #                  color='red', alpha=0.4)
    ax[1,1].plot(range(len(df.ws_avg)), df.ws_avg, color='blue', label='$w^s$')
    # ax[1,1].fill_between(range(len(df.ws_avg)), df.ws_avg + df.ws_std, df.ws_avg - df.ws_std, 
    #                  color='blue', alpha=0.4)
    # ax[1,1].plot(range(len(df.wo_max_avg)), df.wo_max_avg, color='orange', label='$w^o_\max$')
    ax[1,1].set_title('Wage level')
    ax[1,1].legend()

    ax[2,0].plot(range(len(df.dL_avg)), df.dL_avg, color='red', label='all')
    ax[2,0].fill_between(range(len(df.dL_avg)), df.dL_avg + df.dL_std, df.dL_avg - df.dL_std, color='red', alpha=0.4)
    ax[2,0].plot(range(len(df.dL_cp_avg)), df.dL_cp_avg, color='green', label='cp')
    ax[2,0].plot(range(len(df.dL_kp_avg)), df.dL_kp_avg, color='orange', label='kp')
    ax[2,0].set_title('$\Delta L^d$')
    ax[2,0].legend()

    ax[2,1].plot(range(len(df.I_avg)), df.I_avg, color='purple')
    ax[2,1].fill_between(range(len(df.I_avg)), df.I_avg + df.I_std, df.I_avg - df.I_std, color='purple', alpha=0.4)
    ax[2,1].set_title('Income of households')

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

    ax[4,0].plot(range(len(df.EI_avg)), df.EI_avg, label='EI')
    ax[4,0].plot(range(len(df.RS_avg)), df.RS_avg, label='RS')
    ax[4,0].set_title('Investments')
    ax[4,0].legend()

    ax[4,1].plot(range(len(df.avg_pi)), df.avg_pi, label='$\\bar{\pi}$')
    ax[4,1].plot(range(len(df.avg_A)), df.avg_A, label='$\\bar{A}$')
    ax[4,1].plot(range(len(df.avg_B)), df.avg_B, label='$\\bar{B}$')
    ax[4,1].set_title('Productivity')
    ax[4,1].legend()

    ax[5,0].plot(range(len(df.avg_Q_bp)), df.avg_Q_bp, label='bp')
    ax[5,0].plot(range(len(df.avg_Q_lp)), df.avg_Q_lp, label='lp')
    ax[5,0].plot(range(len(df.avg_Q_kp)), df.avg_Q_kp, label='kp')
    ax[5,0].set_title('Production quantity')
    ax[5,0].legend()

    ax[5,1].plot(range(len(df.CPI)), df.CPI)
    ax[5,1].set_title('CPI')

    ax[6,0].plot((range(len(df.bankrupt_bp))), df.bankrupt_bp, label='bp')
    ax[6,0].plot((range(len(df.bankrupt_lp))), df.bankrupt_lp, label='lp')
    ax[6,0].plot((range(len(df.bankrupt_kp))), df.bankrupt_kp, label='kp')
    ax[6,0].legend()
    ax[6,0].set_title('Bankrupty rate')


    plt.tight_layout()
    plt.savefig('plots/macro_ts.png', bbox_inches='tight')


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
    plot_macro_vars()
    plot_income_dist()
    plot_sales_dist()
