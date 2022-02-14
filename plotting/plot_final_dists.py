import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_final_dist():

    df = pd.read_csv('../results/result_data/final_dists.csv')

    # print(df.head())

    fig, [ax1, ax2, ax3]= plt.subplots(1, 3, figsize=(12,5))

    ax1.hist(df.all_I)
    ax1.set_title('Income ($I_{i,t}$)')

    ax2.hist(df.all_w)
    ax2.set_title('Wages ($w_{i,t}$)')

    ax3.hist(df.all_W, bins=100)
    ax3.set_title('Wealth ($W_{i,t}$)')

    plt.tight_layout()

    plt.savefig('plots/final_dist.png')

plot_final_dist()