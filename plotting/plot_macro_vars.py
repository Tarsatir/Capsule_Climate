import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars():

    df = pd.read_csv('../results/result_data/first.csv')

    print(df.head())
    
    fig, ax = plt.subplots(6,2, figsize=(10,15), sharex=True)

    ax[0,0].plot(range(len(df.GDP)), df.GDP, label='total GDP')
    ax[0,0].plot(range(len(df.GDP_I)), df.GDP_I, label='income share')
    ax[0,0].plot(range(len(df.GDP_cp)), df.GDP_cp, label='cp share')
    ax[0,0].plot(range(len(df.GDP_kp)), df.GDP_kp, label='kp share')
    ax[0,0].plot(range(len(df.GDP)), df.Exp_UB, label='UB exp')
    ax[0,0].set_title("GDP")
    ax[0,0].legend()

    ax[0,1].plot(range(len(df.UR)), df.UR)
    ax[0,1].set_title("Unemployment rate")

    ax[1,0].plot(range(len(df.s_avg)), df.s_avg, color='red')
    ax[1,0].fill_between(range(len(df.s_avg)), df.s_avg + df.s_std, df.s_avg - df.s_std, color='red', alpha=0.4)
    ax[1,0].set_title("Savings rate")

    # Wage levels
    ax[1,1].plot(range(len(df.w_avg)), df.w_avg, color='green', label='$\\bar{w}$')
    ax[1,1].fill_between(range(len(df.w_avg)), df.w_avg + df.w_std, df.w_avg - df.w_std, 
                     color='green', alpha=0.4)
    ax[1,1].plot(range(len(df.wr_avg)), df.wr_avg, color='red', label='$w^r$')
    ax[1,1].fill_between(range(len(df.wr_avg)), df.wr_avg + df.wr_std, df.wr_avg - df.wr_std, 
                     color='red', alpha=0.4)
    ax[1,1].plot(range(len(df.ws_avg)), df.ws_avg, color='blue', label='$w^s$')
    ax[1,1].fill_between(range(len(df.ws_avg)), df.ws_avg + df.ws_std, df.ws_avg - df.ws_std, 
                     color='blue', alpha=0.4)
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

    ax[3,1].plot(range(len(df.Deb_tot)), df.Deb_tot, label='total')
    ax[3,1].plot(range(len(df.Deb_cp)), df.Deb_cp, label='cp')
    ax[3,1].plot(range(len(df.Deb_kp)), df.Deb_kp, label='kp')
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


    plt.tight_layout()
    plt.savefig('plots/first.png', bbox_inches='tight')

plot_macro_vars()