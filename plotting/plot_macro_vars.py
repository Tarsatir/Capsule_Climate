import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars():

    df = pd.read_csv('../results/result_data/first.csv')

    print(df.head())
    
    fig, [[ax1, ax2], [ax3, ax4], [ax5, ax6], [ax7, ax8]] = plt.subplots(4,2, figsize=(10,10), sharex=True)

    ax1.plot(range(len(df.GDP)), df.GDP, label='total GDP')
    ax1.plot(range(len(df.GDP_I)), df.GDP_I, label='income share')
    ax1.plot(range(len(df.GDP_cp)), df.GDP_cp, label='cp share')
    ax1.plot(range(len(df.GDP_kp)), df.GDP_kp, label='kp share')
    ax1.plot(range(len(df.GDP)), df.Exp_UB, label='UB exp')
    ax1.set_title("GDP")
    ax1.legend()

    ax2.plot(range(len(df.UR)), df.UR)
    ax2.set_title("Unemployment rate")
    # ax2.fill_between(range(len(df.dL_cp_avg)), df.dL_avg + df.dL_cp_std, df.dL_cp_avg - df.dL_cp_std, color='green', alpha=0.4)

    ax3.plot(range(len(df.s_avg)), df.s_avg, color='red')
    ax3.fill_between(range(len(df.s_avg)), df.s_avg + df.s_std, df.s_avg - df.s_std, color='red', alpha=0.4)
    ax3.set_title("Savings rate")

    ax4.plot(range(len(df.w_avg)), df.w_avg, color='green')
    ax4.fill_between(range(len(df.w_avg)), df.w_avg + df.w_std, df.w_avg - df.w_std, color='green', alpha=0.4)
    ax4.set_title('Wage level')

    ax5.plot(range(len(df.dL_avg)), df.dL_avg, color='red', label='all')
    ax5.fill_between(range(len(df.dL_avg)), df.dL_avg + df.dL_std, df.dL_avg - df.dL_std, color='red', alpha=0.4)
    ax5.plot(range(len(df.dL_cp_avg)), df.dL_cp_avg, color='green', label='cp')
    ax5.plot(range(len(df.dL_kp_avg)), df.dL_kp_avg, color='orange', label='kp')
    ax5.set_title('$\Delta L^d$')
    ax5.legend()

    ax6.plot(range(len(df.I_avg)), df.I_avg, color='purple')
    ax6.fill_between(range(len(df.I_avg)), df.I_avg + df.I_std, df.I_avg - df.I_std, color='purple', alpha=0.4)
    ax6.set_title('Income of households')

    plt.tight_layout()
    plt.savefig('plots/first.png', bbox_inches='tight')

plot_macro_vars()