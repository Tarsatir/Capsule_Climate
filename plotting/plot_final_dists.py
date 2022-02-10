import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_final_dist():

    df = pd.read_csv('../results/result_data/final_dists.csv')

    # print(df.head())

    fig, [ax1, ax2]= plt.subplots(1,2, figsize=(10,5))

    ax1.hist(df.all_I)
    ax1.set_title('Income')

    ax2.hist(df.all_w)
    ax2.set_title('Wages')

    plt.savefig('plots/final_dist.png')

plot_final_dist()