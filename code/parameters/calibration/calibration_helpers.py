import pandas as pd
import numpy as np

def load_data(inputpath, outputpath, n_threads, run_nr=8, return_as_np=True):
    
    # Read input data
    all_dfs = []
    for thread_nr in range(1, n_threads+1):
        df = pd.read_csv(inputpath + f'gsa_input_run{run_nr}_thread{thread_nr}.csv')
        all_dfs.append(df)
    df_input = pd.concat(all_dfs, axis=0, ignore_index=True)

    # Read output data
    all_dfs = []
    for thread_nr in range(1, n_threads+1):
        df = pd.read_csv(outputpath + f'$cal_output_run{run_nr}_thread{thread_nr}.csv')
        all_dfs.append(df)
    df_output = pd.concat(all_dfs, axis=0, ignore_index=True)

    df_input = df_input[df_input["sim_nr"].isin(df_output["sim_nr"])]

    # # Filter only dep vars of interest
    # if dep_vars != None:
    #     df_output = df_output[dep_vars]

    if return_as_np:
        # Convert to numpy array and cut off the simulation number
        # X = df_input.to_numpy()[:, 1:]
        # Y = df_output.to_numpy()

        # run_nrs = np.repeat(df_output['sim_nr'].to_numpy(), 360)

        all_GDP = [col for col in df if col.startswith('GDP')]
        all_U = [col for col in df if col.startswith('U')]
        all_em = [col for col in df if col.startswith('em')]

        gdp = df_output[all_GDP].to_numpy()[:, 1:].flatten()
        # U = df_output[all_U].pct_change(axis=1).fillna(0).to_numpy()[:, 1:].flatten()
        u = df_output[all_U].to_numpy()[:, 1:].flatten()
        # em = df_output[all_em].pct_change(axis=1).fillna(0).to_numpy()[:, 1:].flatten()
        em = df_output[all_em].to_numpy()[:, 1:].flatten()

        # print(df_output[all_U])
        # print(df_output[all_U].pct_change(axis=1).fillna(0).var(axis=1).isna())
        
        # Y = np.array([dGDP, dU, em]).T
        # Y = np.array([gdp, u, em]).T
        # Y = np.array([gdp, em]).T
        Y = np.array([gdp, u]).T

        indepvar = ['κ_upper', 'ω', 'ϵ', 'α_cp', 'p_f', 'prog']
        X = np.repeat(df_input[indepvar].to_numpy(), 360, axis=0)

        return X, Y
    else:
        return df_input, df_output


def handleargs(argv):
    pass