
import numpy as np
from sklearn.covariance import log_likelihood
import statsmodels.api as sm
import scipy

from calibration_helpers import *

class BayesianEstimator:

    def __init__(self, X, Y, Y_true):
        self.X = X
        self.Y = Y
        self.Y_true = Y_true

        self.priors = None
        self.f = None
        self.posterior = None

    def compute_coeff(self):
        
        # Fit likelihood
        self.compute_f()

        self.sample_posterior()


    def compute_f(self):

        # print(self.Y[0,0], self.Y[0,1], self.Y[0,2])
        
        self.f = sm.nonparametric.KDEMultivariateConditional(
                    endog=self.Y,
                    exog=self.X, 
                    dep_type='c' * self.Y.shape[1],
                    indep_type='c' * self.X.shape[1],
                    bw='normal_reference'
                )

    def compute_likelihood(self, x):
        # likelihood = 1
        # for y in self.Y_true:
        #     # print(x)
        #     # print(y)
        #     a = self.f.pdf(endog_predict=y, exog_predict=x)
        #     # print(a)
        #     likelihood *= a
        # return likelihood

        x = np.repeat([x], len(self.Y_true), axis=0)

        return np.product(self.f.pdf(endog_predict=self.Y_true, 
                                     exog_predict=x))

    def checkbounds(self, x_new, bounds):
        for i,x in enumerate(x_new):
            if x < bounds[i][0] or x > bounds[i][1]:
                return True
        return False

    def sample_posterior(self):

        x_curr = np.array([0.05, 0.2, 0.005, 0.8, 0.02, 0.8, 0.5, 0.0])
        bounds = [(0.05, 0.05),
                  (0.2, 0.2),
                  (0.001, 0.01),
                  (0., 1.),
                  (0., 0.05),
                  (0.6, 1.),
                  (0., 1.),
                  (-1., 1.)]
        
        weights = 10 * np.array([0, 0, 0.001, .5, 0.01, .5, .7, .7])
        curr_lh = self.compute_likelihood(x_curr)

        all_x = []

        for _ in range(100):

            selec = np.zeros(len(x_curr))
            selec_idx = np.random.randint(2, len(x_curr))
            selec[selec_idx] = 1.

            epsilon = np.random.normal(0, 0.1, 8)
            x_new = x_curr + epsilon * weights * selec

            while self.checkbounds(x_new, bounds):
                epsilon = np.random.normal(0, 0.1, 8)
                x_new = x_curr + epsilon * weights * selec
            
            new_lh = self.compute_likelihood(x_new)
            print(new_lh)

            alpha = min(1, new_lh / curr_lh)
            print(alpha)

            if np.random.rand() < alpha:
                x_curr = x_new
                curr_lh = new_lh
                all_x.append(x_curr)
                # print(x_curr)

        all_x = np.array(all_x)
        res = {'κ_upper': all_x.T[2], 'ω': all_x.T[3], 'ϵ': all_x.T[4], 
               'α_cp': all_x.T[5], 'p_f': all_x.T[6], 'prog': all_x.T[7]}
        pd.DataFrame(res).to_csv('postdists.csv', index=False)


if __name__ == "__main__":

    inputpath = '../sensitivity/sensitivity_runs/input_data/'
    outputpath = 'calibrationdata/'
    n_threads = 4

    # Y_true = {
    #     'GDP_1st': 0.482,
    #     'GDP_2nd': 0.47,
    #     'U_1st': 0.06097,
    #     'em2030': 112.72,
    #     'em2030': 121.82,
    #     'em2040': 127.27
    # }

    # dep_vars = Y_true.keys()

    Y_true = pd.read_csv(outputpath + 'truedata.csv')[["dGDP", "em"]].to_numpy()

    X, Y = load_data(inputpath, outputpath, n_threads)

    estimator = BayesianEstimator(X, Y, Y_true)
    estimator.compute_coeff()
