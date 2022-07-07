import pandas as pd

def load_data(inputpath, outputpath, n_threads, return_as_np=True, dep_vars=None):
    """_summary_

    Args:
        inputpath (String): _description_
        outputpath (String): _description_
        n_threads (Int): _description_

    Returns:
        _type_: _description_
    """
    
    # Read input data
    all_dfs = []
    for thread_nr in range(1, n_threads+1):
        df = pd.read_csv(inputpath + f'gsa_input_run7_thread{thread_nr}.csv')
        all_dfs.append(df)
    df_input = pd.concat(all_dfs, axis=0, ignore_index=True)

    # Read output data
    all_dfs = []
    for thread_nr in range(1, n_threads+1):
        df = pd.read_csv(outputpath + f'gsa_output_run7_thread{thread_nr}.csv')
        all_dfs.append(df)
    df_output = pd.concat(all_dfs, axis=0, ignore_index=True)

    df_input = df_input[df_input["sim_nr"].isin(df_output["sim_nr"])]

    # Filter only dep vars of interest
    if dep_vars != None:
        df_output = df_output[dep_vars]

    if return_as_np:
        # Convert to numpy array and cut off the simulation number
        X = df_input.to_numpy()[:, 1:]
        Y = df_output.to_numpy()

        return X, Y
    else:
        return df_input, df_output


def handleargs(argv):
    pass