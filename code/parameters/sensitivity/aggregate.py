
import pandas as pd
import numpy as np
import os 


def aggregate_data(data_path):
    # Read in the data
    df = pd.read_csv(data_path)
    # Group by 'sim_nr' and then apply a lambda function to take the mean of the last 20 rows for each column
    grouped_df = df.groupby('sim_nr').apply(lambda x: x[-20:].mean(numeric_only=True))
    
    # Reset the index while avoiding an extra 'sim_nr' column
    grouped_df.reset_index(drop=True, inplace=True)
    
    # Add the 'sim_nr' values as a new column at the beginning of the DataFrame
    grouped_df['sim_nr'] = df['sim_nr'].unique()
    
    # Reorder the columns to make 'sim_nr' the first column
    cols = ['sim_nr'] + [col for col in grouped_df.columns if col != 'sim_nr']
    grouped_df = grouped_df[cols]

    grouped_df.to_csv(f'sensitivity_runs/output_data_agg/{data_path.split("/")[-1]}', index=False)  
    return grouped_df



# Converting the Julia function to a Python function
def get_output_path(parl_id: int, run_nr: int) -> str:
    return f"sensitivity_runs/output_data/gsa_output_run{run_nr}_thread{parl_id}.csv"


X_labels = {
    "α_maxdev": [0.005, 0.5],
    "ρ": [0.05, 0.8],
    "prog": [-1.0, 1.0],
    "μ1": [0.0, 0.5],
    "ω": [0.0, 1.0],
    "λ": [0.0, 1.0],
    "ϵ_w": [0.0, 0.1],
    "ϵ_μ": [0.0, 0.1],
    "κ_upper": [0.0, 0.01],
    "ψ_E": [0.0, 0.25],
    "ψ_Q": [0.0, 0.25],
    "ψ_P": [0.0, 0.25],
    "p_f": [0.0, 1.0]
}
run_nr = 9






if __name__ == "__main__":
    csv_files = [f for f in os.listdir("sensitivity_runs/output_data") if f.endswith('.csv')]

    # Loop over each csv file
    for idx, csv_file in enumerate(csv_files):
        # Use the index of the csv file as its number
        csv_filenumber = idx + 1  # Adding 1 because index starts from 0
        
        out_path = get_output_path(csv_filenumber, run_nr)
        new_df = aggregate_data(out_path)
    