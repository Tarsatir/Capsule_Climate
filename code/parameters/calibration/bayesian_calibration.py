
import numpy as np
# from sklearn.covariance import log_likelihood
import statsmodels.api as sm
# import scipy
import matplotlib.pyplot as plt
import julia

# import pymc3 as pm


from calibration_helpers import *

class BayesianEstimator:

    def __init__(self, df_X, df_Y, y_true):
        self.df_X = df_X
        self.df_Y = df_Y
        self.y_true = y_true
        self.n_samples = len(df_X)

        self.get_col_labels()

    def get_col_labels(self):

        self.all_gdp = [col for col in self.df_Y if col.startswith('GDP')]
        self.all_u = [col for col in self.df_Y if col.startswith('U')]
        self.all_em = [col for col in self.df_Y if col.startswith('em')]

    def build_Julia_session(self):
        jl = julia.Julia(compiled_modules=False)
        jl.include("bayesian_calibration.jl")
        return jl

    def sample_coeff(self):

        # Make a Julia session
        # jl = self.build_Julia_session()
        
        # # Initialize with first observed combination
        # x_curr = self.get_x(0)
        # # y = self.get_y(0)
        
        # gdp, u, em = jl.run_model(x_curr[1:])
        # y = np.array([gdp, u, em]).T
        # curr_lh = self.compute_likelihood(y)

        # # print(curr_lh)

        # # min_x = x_curr
        # # min_y = y
        # # min_lh = curr_lh

        # for i_row in range(1, self.n_samples):

        #     # Gather data and transform
        #     y = self.get_y(i_row)

        #     # print(self.y_true)
        #     # print(y)
        #     # return

        #     # Compute likelihood of observed series y
        #     new_lh = self.compute_likelihood(y)
        #     print(i_row, new_lh)

        #     # alpha = np.min([1., new_lh / curr_lh])
        #     # if np.random.rand() < alpha:
        #     #     x_curr = self.get_x(i_row)
        #     #     curr_lh = new_lh
        #     #     print(f'accepted at alpha = {alpha}')

        #     if new_lh > min_lh:
        #         min_x = self.get_x(i_row)
        #         min_y = y
        #         min_lh = new_lh

        # print(min_x)
        # print(min_y)
        # print(min_lh)

        # fig, ax = plt.subplots(1,3)
        # ax[0].hist(self.y_true[:, 0], alpha=0.5, density=True, bins=20)
        # ax[0].hist(min_y[:, 0], alpha=0.5, density=True, bins=20)

        # ax[1].hist(self.y_true[:,1], alpha=0.5, density=True, bins=20)
        # ax[1].hist(min_y[:,1], alpha=0.5, density=True, bins=20)

        # ax[2].hist(self.y_true[:,2], alpha=0.5, density=True, bins=20)
        # ax[2].hist(min_y[:,2], alpha=0.5, density=True, bins=20)

        # plt.show()

    def get_y(self, i):

        df_y = self.df_Y.iloc[[i]]

        gdp = df_y[self.all_gdp].to_numpy()[:, 1:].flatten()
        u = df_y[self.all_u].to_numpy()[:, 1:].flatten() * 100
        em = df_y[self.all_em].to_numpy()[:, 1:].flatten()

        return np.array([gdp, u, em]).T

    def get_x(self, i):
        return self.df_X.iloc[[i]].to_numpy()[0]

    def compute_likelihood(self, y):

        # Fit function f
        f = sm.nonparametric.KDEMultivariate(
                data=y, 
                var_type='c'*y.shape[1],
                bw='normal_reference'
            )
        likelihood = np.product(f.pdf(self.y_true))
        return likelihood


if __name__ == "__main__":

    inputpath = '../sensitivity/sensitivity_runs/input_data/'
    outputpath = 'calibrationdata/'
    n_threads = 16
    run_nr = 9

    Y_true = pd.read_csv(outputpath + 'truedata.csv')[["dGDP", "Ur", "em"]].to_numpy()

    df_X, df_Y = load_data(inputpath, outputpath, n_threads, run_nr=run_nr, return_as_np=False)

    estimator = BayesianEstimator(df_X, df_Y, Y_true)
    estimator.sample_coeff()
