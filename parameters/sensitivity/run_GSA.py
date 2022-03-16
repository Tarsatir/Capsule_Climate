"""
This file contains code to run PAWN global sensitivity analysis.

"""

from SAFEpython import PAWN
import scipy.stats as stats
import numpy as np

from SAFEpython.sampling import AAT_sampling
import SAFEpython.plot_functions as pf # module to visualize the results

import matplotlib.pyplot as plt


def call_AAT_sampling(samp_strat, M, X_labels, N):
    """
    Samples parameters to use for simulation runs.
    """
    
    # Define distribution of parameters
    distr_fun = [stats.uniform] * M

    distr_par = [np.nan] * M
    for i,key in enumerate(X_labels):
        distr_par[i] = [X_labels[key][0], X_labels[key][1] - X_labels[key][0]]

    X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
    
    return X


def run_PAWN(X_labels, X, Y, n=10):
    """
    Runs code required for PAWN sensitivity analysis.
    """

    KS_median, KS_mean, KS_max = PAWN.pawn_indices(X, Y, n)

    plt.figure()
    pf.boxplot1(KS_max, X_Labels=X_labels, Y_Label='Ks (max)')
    plt.savefig('parameters/sensitivity/sensitivity_runs/sensitivity_plot1.png')