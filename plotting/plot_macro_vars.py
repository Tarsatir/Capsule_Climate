import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars():

    df = pd.read_csv('../results/result_data/first.csv')

    print(df.head())
    
    fig, [ax1, ax2, ax3, ax4] = plt.subplots(4,1, figsize=(7,10), sharex=True)

    ax1.plot(range(len(df.GDP)), df.GDP, label='total GDP')
    ax1.plot(range(len(df.GDP_I)), df.GDP_I, label='income share')
    ax1.plot(range(len(df.GDP_cp)), df.GDP_cp, label='cp share')
    ax1.plot(range(len(df.GDP_kp)), df.GDP_kp, label='kp share')
    ax1.set_title("GDP")
    ax1.legend()

    ax2.plot(range(len(df.UR)), df.UR)
    ax2.set_title("Unemployment rate")

    ax3.plot(range(len(df.s_avg)), df.s_avg, color='red')
    ax3.fill_between(range(len(df.s_avg)), df.s_avg + df.s_std, df.s_avg - df.s_std, color='red', alpha=0.4)
    ax3.set_title("Savings rate")

    ax4.plot(range(len(df.w_avg)), df.w_avg, color='green')
    ax4.fill_between(range(len(df.w_avg)), df.w_avg + df.w_std, df.w_avg - df.w_std, color='green', alpha=0.4)
    ax4.set_title('Wage level')

    plt.tight_layout()
    plt.savefig('plots/first.png', bbox_inches='tight')

plot_macro_vars()