"""
This file contains code to run PAWN global sensitivity analysis.

"""

from SAFEpython import PAWN
import scipy.stats as stats
import numpy as np

from SAFEpython.sampling import AAT_sampling
import SAFEpython.plot_functions as pf # module to visualize the results
from SAFEpython.util import aggregate_boot  # function to aggregate the bootstrap results

import matplotlib.pyplot as plt
from matplotlib import rc


def call_AAT_sampling(samp_strat, X_labels, N):
    """
    Samples parameters to use for simulation runs.
    """
    
    # Define distribution of parameters
    distr_fun = [stats.uniform] * len(X_labels)

    distr_par = [np.nan] * len(X_labels)
    for i,key in enumerate(X_labels):
        distr_par[i] = [X_labels[key][0], X_labels[key][1] - X_labels[key][0]]

    X = AAT_sampling(samp_strat, len(X_labels), distr_fun, distr_par, N)
    
    return X


def run_PAWN(X_labels, X, Y, type, run_nr, name_dep_var, n=10, Nboot=3000):
    """
    Runs code required for PAWN sensitivity analysis.
    """

    KS_median, KS_mean, KS_max = PAWN.pawn_indices(X, Y, n, Nboot=Nboot)

    # KS_median_m, KS_median_lb, KS_median_ub = aggregate_boot(KS_median) # shape (M,)
    KS_mean_m, KS_mean_lb, KS_mean_ub = aggregate_boot(KS_mean) # shape (M,)
    # KS_max_m, KS_max_lb, KS_max_ub = aggregate_boot(KS_max) # shape (M,)

    X_labels = [f'${l}$' for l in X_labels]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # Plot bootstrapping results (for instance for KS_max):
    # plt.figure()
    # pf.boxplot1(KS_max_m, S_lb=KS_max_lb, S_ub=KS_max_ub,
    #             X_Labels=X_labels, Y_Label=f'{name_dep_var} (max)')
    # plt.savefig(f'parameters/sensitivity/sensitivity_runs/sa_{type}_max_plot{run_nr}.png')

    plt.figure(figsize=(8, 3))
    pf.boxplot1(KS_mean_m, S_lb=KS_mean_lb, S_ub=KS_mean_ub,
                X_Labels=X_labels, Y_Label=f'{name_dep_var}')
    plt.tight_layout()
    plt.savefig(f'parameters/sensitivity/sensitivity_plots/sa_{type}_mean_plot{run_nr}.pdf')

    # YF, FU, FC, xc = PAWN.pawn_plot_cdf(X, Y, n, cbar=True, n_col=3, labelinput=X_labels)
    # plt.show()

    # KS = PAWN.pawn_plot_ks(YF, FU, FC, xc)
    # plt.show()

    # plt.figure()
    # pf.boxplot1(KS_median_m, S_lb=KS_median_lb, S_ub=KS_median_ub,
    #             X_Labels=X_labels, Y_Label=f'{name_dep_var} (median)')
    # plt.savefig(f'parameters/sensitivity/sensitivity_runs/sensitivity_median_plot{run_nr}.png')