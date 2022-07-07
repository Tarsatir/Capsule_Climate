
import numpy as np
import statsmodels.api as sm
import scipy

from calibration_helpers import *

class BayesianEstimator:

    def __init__(self, X, Y, Y_true):
        self.X = X
        self.Y = Y
        self.Y_true = list(Y_true.values())

        self.priors = None
        self.likelihood = None
        self.posterior = None

    def compute_coeff(self):
        
        # Fit likelihood
        self.compute_likelihood()

    def compute_likelihood(self):
        
        self.likelihood = sm.nonparametric.KDEMultivariateConditional(
                            endog=self.Y,
                            exog=self.X, 
                            dep_type='c' * self.Y.shape[1],
                            indep_type='c' * self.X.shape[1],
                            bw='normal_reference'
                        )


if __name__ == "__main__":

    inputpath = '../sensitivity/sensitivity_runs/input_data/'
    outputpath = '../sensitivity/sensitivity_runs/output_data/'
    n_threads = 16

    Y_true = {
        'GDP_1st': 0.482,
        'GDP_2nd': 0.47,
        'U_1st': 0.06097,
        'em2030': 112.72,
        'em2030': 121.82,
        'em2040': 127.27
    }

    dep_vars = Y_true.keys()

    X, Y = load_data(inputpath, outputpath, n_threads, dep_vars=dep_vars)

    estimator = BayesianEstimator(X, Y, Y_true)
    estimator.compute_coeff()
