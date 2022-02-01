import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_macro_vars():

    df = pd.read_csv('../results/result_data/first.csv')
    
    fig, [ax1, ax2] = plt.subplots(2,1, figsize=(7,7), sharex=True)

    ax1.plot(range(len(df.GDP)), df.GDP)
    ax1.set_title("GDP")

    ax2.plot(range(len(df.UR)), df.UR)
    ax2.set_title("Unemployment rate")

    plt.tight_layout()
    plt.savefig('plots/first.png')

plot_macro_vars()