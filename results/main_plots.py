


from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
import numpy as np
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
#import get_cmap
from matplotlib.cm import get_cmap
from pybdm import BDM
import math



cmap = get_cmap('viridis')





def Create_dataframe(path1):
    warnings.simplefilter(action='ignore', category=FutureWarning)

    path = path1  
    #print files in path
    print(os.listdir(path))
    all_files = glob.glob(os.path.join(path, "*.csv"))
    all_files.sort(key=lambda x: float(x.split('_')[1]))
    start_step = float(all_files[0].split('_')[1])
    end_step = float(all_files[-1].split('_')[1])
    middle_number = float(all_files[0].split('_')[1])  # Assuming the list is already sorted
    runs = sum(1 for file in all_files if float(file.split('_')[1]) == middle_number)

    #create linespace 
    steps= int(len(all_files)/runs)
    #print(steps)
    linspace = np.linspace(float(start_step), float(end_step), steps)


    df_from_each_file = (pd.read_csv(f) for f in all_files)
    df = pd.concat(df_from_each_file, ignore_index=True)

    df_list = np.array_split(df, len(all_files))

    df_mean = pd.DataFrame()
    df_std = pd.DataFrame()
    for i in range(0, len(all_files)):
        # df_mean = df_mean.append(pd.DataFrame(df_list[i].iloc[-10:].mean()).transpose())
        # df_std = df_std.append(pd.DataFrame(df_list[i].iloc[-10:].std()).transpose()) 
        mean_values = pd.DataFrame(df_list[i].iloc[-10:].mean()).transpose()
        std_values = pd.DataFrame(df_list[i].iloc[-10:].std()).transpose()
        
        # Use pd.concat() instead of append
        df_mean = pd.concat([df_mean, mean_values], ignore_index=True)
        df_std = pd.concat([df_std, std_values], ignore_index=True)

    df_mean_mean = pd.DataFrame()
    df_mean_std = pd.DataFrame()
    
    # Create lists to collect DataFrames before concatenation
    mean_mean_list = []
    mean_std_list = []

    for i in range(0, steps):
        mean_values = pd.DataFrame(df_mean.iloc[runs * i : runs * i + runs].mean()).transpose()
        std_values = pd.DataFrame(df_mean.iloc[runs * i : runs * i + runs].std()).transpose()
        
        # Add to the respective lists
        mean_mean_list.append(mean_values)
        mean_std_list.append(std_values)

    # Concatenate the collected DataFrames
    df_mean_mean = pd.concat(mean_mean_list, ignore_index=True)
    df_mean_std = pd.concat(mean_std_list, ignore_index=True)
    if len(linspace) != runs:
        print("Warning: len(linspace) is not equal to runs")
        print("len(linspace) = ", len(linspace))
        print("runs = ", runs)
    
    df_mean_mean.index = linspace[::-1]
    df_mean_std.index = linspace[::-1]

    return df_mean_mean, df_mean_std, linspace

def prog_exp(ax):
    cmap = get_cmap('viridis')
    path1 = "/22 Data"
    path2 = "/23 Data"
    # hightax = os.getcwd() + path1
    # lowtax = os.getcwd() + path2
    carb_tax = 0.5
    no_tax = 0.1



    df_mean_mean1, df_mean_std1, linspace1 = Create_dataframe(path1) #with carbon tax
    df_mean_mean2, df_mean_std2, linspace2 = Create_dataframe(path2)
    #normalize gdp by first entry
    #show me the head of df_mean_mean1
    #print(df_mean_mean1.head(10))
    #check if "GDP" exists in df_mean_mean1

    scaling_factor = df_mean_mean1['GDP'].iloc[-1] #with carbon tax
    scaling_factor2 = df_mean_mean2['GDP'].iloc[-1]

    df_mean_mean1['GDP'] = df_mean_mean1['GDP']/scaling_factor#with carbon tax
    df_mean_std1['GDP'] = df_mean_std1['GDP']/scaling_factor #with carbon tax
    df_mean_mean2['GDP'] = df_mean_mean2['GDP']/scaling_factor2
    df_mean_std2['GDP'] = df_mean_std2['GDP']/scaling_factor2


    dn1 = 'GDP'

    ax.errorbar(df_mean_mean1.index, df_mean_mean1[dn1], yerr=df_mean_std1[dn1], capsize=2, fmt='s', markersize=2, color=cmap(carb_tax), elinewidth=0.5) #with carbon tax
    ax.set_xlabel('Progressivity')
    ax.set_ylabel('Normalized GDP \n (with carbon tax)')

    # Create a twin Axes sharing the x-axis
    ax2 = ax.twinx()
    ax2.errorbar(df_mean_mean2.index, df_mean_mean2[dn1], yerr=df_mean_std2[dn1], capsize=2, fmt='s', markersize=2, color=cmap(no_tax), elinewidth=0.5)
    ax2.set_ylabel('Normalized GDP \n(without carbon tax)')
    

    # Adjust the scale of the second y-axis
    #ax2.set_yscale('linear')

    # legend_handles = [
    #     mpatches.Patch(color=cmap(carb_tax), label='0.0 Carbon Tax'),
    #     mpatches.Patch(color=cmap(no_tax), label='0.5 Carbon Tax')
    # ]
    ax.yaxis.label.set_color(cmap(carb_tax))
    ax2.yaxis.label.set_color(cmap(no_tax))
    #adjust ax2 y axis limit to match ax1
    ax2.set_ylim(ax.get_ylim())
    ax2.tick_params(axis='y', which='both', labelright=False, right=False)
    
    #axins = inset_axes(ax2, width="30%", height="30%", loc=2)
    axins = inset_axes(ax2, width="30%", height="30%", loc='upper left', bbox_to_anchor=(0.1, 0, 1, 1), bbox_transform=ax2.transAxes)
    axins.errorbar(df_mean_mean2.index, df_mean_mean2[dn1], yerr=df_mean_std2[dn1], capsize=2, fmt='s', markersize=2, color=cmap(no_tax), elinewidth=0.5)
    axins.set_adjustable('datalim')
    axins.autoscale_view()
    # Add the custom legend to the figure
    #plt.legend(handles=legend_handles, loc='upper left', prop={'size': 7})


def plot_GDP_emissions2(ax, df_mean, df_std, linspace, position):
    cmap = get_cmap('viridis')
    em = 0.75
    gdp = 0.5
    dn1 = 'GDP'
    dn2 = 'em_index'
    
    if position == 'left':
        ax.errorbar(linspace, df_mean[dn1], yerr=df_std[dn1], capsize=2, fmt='s', markersize=2, color=cmap(gdp), elinewidth=0.5, label=dn1)
        ax.set_ylabel('GDP')
        ax.set_xlabel('Carbon Tax \n (High Green Capacity)')
        #show legend 
        #ax.legend(dn1, loc='upper center')
        #ax.tick_params(axis='x', which='both', labelbottom=False, bottom=False)
        #ax.yaxis.get_offset_text().set_fontsize(3)
        
        ax2 = ax.twinx()
        ax2.errorbar(linspace, df_mean[dn2], yerr=df_std[dn2], capsize=2, fmt='o', markersize=2, color=cmap(em), elinewidth=0.5, label=dn2)
        #ax2.tick_params(axis='y', which='both', labelright=False, right=False)
        #rename legent and show it


        

    elif position == 'right':
        ax.errorbar(linspace, df_mean[dn1], yerr=df_std[dn1], capsize=2, fmt='s', markersize=2, color=cmap(gdp), elinewidth=0.5, label=dn1)
        #ax.tick_params(axis='y', which='both', labelleft=False, left=False)
        ax.set_xlabel('Carbon Tax \n (Low Green Capacity)')
        #ax.tick_params(axis='x', which='both', labelbottom=False, bottom=False)
        
        ax2 = ax.twinx()
        ax2.errorbar(linspace, df_mean[dn2], yerr=df_std[dn2], capsize=2, fmt='o', markersize=2, color=cmap(em), elinewidth=0.5, label=dn2)
        ax2.set_ylabel('Final CO2 Emission-Index')

        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines + lines2, labels + labels2, loc='lower center')



def plot_new_panel():

    path1 = "../data/OFAT/"
    df_mean_mean1, df_mean_std1, linspace1 = Create_dataframe(path1)

    fig = plt.figure(figsize=(4, 4))  # Adjust size as needed
    spec = GridSpec(ncols=1, nrows=1, figure=fig)

    # Upper left plot
    ax1 = fig.add_subplot(spec[0, 0])

    ax1.text(-0.1, 1.1, 'a)', transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top')

    # Use science plots
    #plt.style.use(['science', 'ieee'])
    # Call your plot function or other plotting code here
    plot_GDP_emissions2(ax1, df_mean_mean1, df_mean_std1, linspace1, position='left')

    # Make background white
    fig.patch.set_facecolor('white')
    plt.savefig('../results/critical_transition_a.pdf', dpi=300)
    #plt.show()






def create_file_list(folder_path):
    # Get a list of all files in the folder
    all_files = os.listdir(folder_path)

    # Filter for only .csv files
    csv_files = [file for file in all_files if file.endswith('.csv')]

    # Prepend the folder path to each filename
    full_paths = [os.path.join(folder_path, file) for file in csv_files]

    return full_paths

def compute_bdm(df, variables, cut=20):

    def binarize(x):
        return np.where(x < np.mean(x), 0, 1)

    def func_bdm(x):
        bdm = BDM(ndim=1)
        return bdm.bdm(x, normalized=True)

    def dataprep(x, cut=20):
        dat = np.array(x)
        dat = np.diff(dat)
        dat = [list(dat[i:i+cut]) for i in range(0, len(dat), cut)]
        dat = [binarize(x) for x in dat]
        return dat

    results = pd.DataFrame()
    for var in variables:
        dat = dataprep(df[var], cut)
        KC = [func_bdm(x) for x in dat]
        KC = pd.Series(KC).rolling(10).mean()
        results[var] = KC

    results.index = results.index * 20

    return results

def compute_averages_and_std(file_list,variables):
    # Read the first file to get the number of rows and column names
    initial_df = pd.read_csv(file_list[0])
    n_rows = len(initial_df)
    column_names = initial_df.columns

    # Initialize a dataframe to store sums, squares, and counts
    sum_df = pd.DataFrame(0.0, index=range(n_rows), columns=column_names)
    squares_df = pd.DataFrame(0.0, index=range(n_rows), columns=column_names)
    count_df = pd.DataFrame(0, index=range(n_rows), columns=column_names)

    # Iterate over each file
    for filename in file_list:
        # Load the file into a dataframe
        df = pd.read_csv(filename)
        #print dimension of df
        #print(df.shape)
        # Compute the KC values for this file
        KC_df = compute_bdm(df,variables)
        #print(KC_df.shape)
        # Add the values to the sum dataframe and squares dataframe
        sum_df += KC_df
        squares_df += KC_df**2

        # Update counts (this handles potential differences in timesteps between files)
        count_df += df.notnull()

    # Compute averages and standard deviations
    avg_df = sum_df / count_df
    std_df = np.sqrt((squares_df / count_df) - (avg_df ** 2))
    #remove columns with only NaN values
    avg_df = avg_df.dropna(axis=1, how='all')
    std_df = std_df.dropna(axis=1, how='all')

    return avg_df, std_df



def plot_on_ax(ax1, file_list, variable_names):
    mean, std = compute_averages_and_std(file_list, variable_names)
    mean = mean.dropna(axis=0, how='all')
    std = std.dropna(axis=0, how='all')
    
    palette = sns.color_palette("mako_r", len(variable_names))
    
    # Plot the dataframe on ax1
    mean.plot(ax=ax1, color=palette, linewidth=1.0)
    
    # Plot std on ax1 for all variables
    for i in range(len(variable_names)):
        ax1.fill_between(std.index, mean[variable_names[i]] - std[variable_names[i]], mean[variable_names[i]] + std[variable_names[i]], alpha=0.2, color=palette[i])
    
    ax2 = ax1.twinx()
    df_carbon = pd.read_csv(file_list[0])
    ax2.plot(df_carbon['τᶜ_ts'][200:], color='grey', linestyle='--')
    ax2.set_ylabel('Carbon Tax')

    ax1.set_xlabel('Time')
    ax1.set_ylabel('Normalized Kolmogorov Complexity')
    
    ax1.legend(loc='lower right', fontsize='small')

def compute_simple_averages(file_list, variables):
    # Initialize sum_df and count_df with zeros based on the shape of the first file
    initial_df = pd.read_csv(file_list[0])
    n_rows = len(initial_df)
    sum_df = pd.DataFrame(0.0, index=range(n_rows), columns=variables)
    count_df = pd.DataFrame(0, index=range(n_rows), columns=variables)

    # Iterate over each file
    for filename in file_list:
        # Load the file into a dataframe
        df = pd.read_csv(filename)

        # Only take the columns specified by 'variables'
        df_filtered = df[variables]

        # Add the values to sum_df and update count_df
        sum_df += df_filtered
        count_df += df_filtered.notnull()

    # Compute averages
    avg_df = sum_df / count_df
    
    # Remove columns with only NaN values
    avg_df = avg_df.dropna(axis=1, how='all')
    
    return avg_df

def calculate_susceptibility(df, time_window):
    # Calculate the rate of change for avg_norm_gdp and avg_norm_emission
    df['delta_gdp'] = df['GDP'].diff()
    df['delta_emission'] = df['em_index'].diff()

    # Drop NA rows caused by differencing
    df = df.dropna()

    # Initialize an empty DataFrame to store susceptibility values
    susceptibility_df = pd.DataFrame()
    print(len(df))
    # Loop through DataFrame based on time_window
    for start in range(0, len(df) - time_window + 1):
        end = start + time_window

        # Create a slice of the DataFrame for the current time window
        df_slice = df.iloc[start:end]

        # Calculate susceptibility for this time window
        susceptibility = (df_slice['delta_emission'] / df_slice['delta_gdp']).mean()

        # Add susceptibility value to the output DataFrame
        susceptibility_df.loc[df_slice.index[-1], 'susceptibility'] = susceptibility
    print(len(susceptibility_df))

    return susceptibility_df

def plot_susceptibility(ax, file_list, variable_names, time_window=5):
    mean = compute_simple_averages(file_list, variable_names)
    mean = (mean - mean.mean()) / mean.std()
    result = calculate_susceptibility(mean, time_window)
 
    ax.plot(result.index[200:], result['susceptibility'][200:], label='Susceptibility')
    ax.set_ylabel('Susceptibility')
    ax.legend()
    
    
    return ax

def plot_KC(data_path='../data/multirun/', save_path='../results/KC_plot.pdf'):
    """
    Creates and saves a KC plot with specified parameters.

    Parameters:
        data_path (str): Path to the directory containing the data files.
        save_path (str): Path where the plot will be saved.

    Returns:
        None
    """
    # Generate the file list
    file_list = create_file_list(data_path)

    # Define the variable names
    variable_names1 = ['GDP', 'em_index']
    variable_names2 = ['U', 'em_index', 'debt_kp']

    # Create the figure object
    fig = plt.figure(figsize=(6, 4))

    # Define GridSpec layout
    gs = GridSpec(2, 1, height_ratios=[1, 0.3])

    # Create the subplots (ax1 is the larger plot, ax2 is the smaller one)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)  # sharex ensures the x-axis is aligned

    # Hide x-axis label for the top plot
    ax1.set_xlabel('')
    ax1.xaxis.set_tick_params(labelbottom=False)

    # Plot data on the axes
    plot_on_ax(ax1, file_list, variable_names2)
    plot_susceptibility(ax2, file_list, variable_names1)

    # Adjust layout
    plt.tight_layout()

    # Add figure-level text and formatting
    fig.text(0.5, 0, 'Susceptibility over Time', ha='center')
    ax1.set_xlim(left=200, right=ax1.get_xlim()[1])
    ax2.set_xlim(left=200, right=ax2.get_xlim()[1])

    # Add vertical line at x = 664
    ax1.axvline(x=664, color='black', linestyle='--', linewidth=0.2)
    ax2.axvline(x=664, color='black', linestyle='--', linewidth=0.2)

    # Add shaded region between x = 200 and x = 360
    ax1.axvspan(xmin=200, xmax=360, color='gray', alpha=0.2)
    ax2.axvspan(xmin=200, xmax=360, color='gray', alpha=0.2)

    # Add white background
    fig.patch.set_facecolor('white')

    # Save the figure
    plt.savefig(save_path, dpi=300)
    plt.close(fig)



def process_household_data(household_data_path, cpi_data_path, num_slices=3):
    """
    Process household data and CPI data to compute average income and wealth for percentiles.

    Parameters:
    - household_data_path: str, path to the household data CSV file.
    - cpi_data_path: str, path to the CPI data CSV file.
    - num_slices: int, number of percentiles to divide households into (default: 3).

    Returns:
    - df_avg_all_I: pd.DataFrame, DataFrame containing average income data for percentiles.
    - df_avg_all_W: pd.DataFrame, DataFrame containing average wealth data for percentiles.
    """
    # Load data
    df_combined = pd.read_csv(household_data_path)
    df_CPI = pd.read_csv(cpi_data_path)

    # Extract base file name for saving results
    household_filename = os.path.basename(household_data_path).split('.')[0]
    
    # Assume the second column of df_CPI contains the necessary CPI values
    df_CPI = df_CPI.iloc[:, 1].reset_index(drop=True)

    # Sort the data for initial percentile assignment
    df0 = df_combined[df_combined['timestamp'] == 1]  # Use timestamp 1 for initial percentile calculations
    df0_I = df0.sort_values(by=['all_I'])
    df0_w = df0.sort_values(by=['all_W'])

    # Save household IDs of each percentile
    id_I = df0_I['hh_id'].tolist()
    id_w = df0_w['hh_id'].tolist()
    percentiles_I = []
    percentiles_W = []

    slice_size = len(id_I) // num_slices

    for i in range(num_slices):
        start_idx = i * slice_size
        end_idx = (i + 1) * slice_size
        percentiles_I.append(id_I[start_idx:end_idx])
        percentiles_W.append(id_w[start_idx:end_idx])

    # Initialize dictionaries to store the average of all_I and all_W for each percentile
    avg_all_I = {f'Percentile_{i+1}': [] for i in range(num_slices)}
    avg_all_W = {f'Percentile_{i+1}': [] for i in range(num_slices)}

    # Process each timestamp in the combined CSV file
    for timestamp in df_combined['timestamp'].unique():
        df = df_combined[df_combined['timestamp'] == timestamp]

        # Ensure the timestamp matches a valid CPI index
        cpi_index = timestamp - 1  # Convert 1-based timestamp to 0-based index
        if cpi_index < 0 or cpi_index >= len(df_CPI):
            print(f"Skipping timestamp {timestamp}: Out of bounds for CPI data")
            continue

        for i, (p_I, p_W) in enumerate(zip(percentiles_I, percentiles_W)):
            avg_I = df[df['hh_id'].isin(p_I)]['all_I'].mean()
            avg_W = df[df['hh_id'].isin(p_W)]['all_W'].mean()

            avg_I = avg_I / df_CPI.iloc[cpi_index]
            avg_W = avg_W / df_CPI.iloc[cpi_index]

            percentile_key = f'Percentile_{i+1}'
            avg_all_I[percentile_key].append(avg_I)
            avg_all_W[percentile_key].append(avg_W)

    # Convert the dictionaries to DataFrames
    df_avg_all_I = pd.DataFrame.from_dict(avg_all_I, orient='index').transpose()
    df_avg_all_W = pd.DataFrame.from_dict(avg_all_W, orient='index').transpose()

    # Save the DataFrames to CSV files
    output_folder = os.path.dirname(household_data_path)
    avg_all_I_path = os.path.join(output_folder, f"{household_filename}_avg_I.csv")
    avg_all_W_path = os.path.join(output_folder, f"{household_filename}_avg_W.csv")
    
    df_avg_all_I.to_csv(avg_all_I_path, index=False)
    df_avg_all_W.to_csv(avg_all_W_path, index=False)

    print(f"Files saved:\n- {avg_all_I_path}\n- {avg_all_W_path}")

    return df_avg_all_I, df_avg_all_W



def generate_incomeshock_plot(
    no_shock_income_path, no_shock_wealth_path,
    up_shock_income_path, up_shock_wealth_path,
    down_shock_income_path, down_shock_wealth_path,
    output_path, smooth=True, window_size=10, dpi_value=600
):
    """
    Generate and save a plot for relative income and wealth changes due to shocks.

    Parameters:
    - no_shock_income_path: str, path to the no-shock income CSV file.
    - no_shock_wealth_path: str, path to the no-shock wealth CSV file.
    - up_shock_income_path: str, path to the +50% shock income CSV file.
    - up_shock_wealth_path: str, path to the +50% shock wealth CSV file.
    - down_shock_income_path: str, path to the -50% shock income CSV file.
    - down_shock_wealth_path: str, path to the -50% shock wealth CSV file.
    - output_path: str, where the resulting plot should be saved.
    - smooth: bool, whether to apply rolling smoothing (default: True).
    - window_size: int, window size for rolling mean smoothing (default: 10).
    - dpi_value: int, resolution of the saved plot (default: 600).
    """
    vmap = get_cmap('viridis')

    def get_percentile_group(col_name):
        # Extract the last string (assuming it's a number) from the column name
        percentile_str = col_name.split("_")[-1]  # Adjust based on your column naming convention
        try:
            percentile = int(percentile_str)
        except ValueError:
            # Handle cases where the conversion to int might fail
            return "Unknown Group"
        
        # Assigning labels based on the percentile value
        if percentile <= 1:
            return "Lowest "
        elif percentile <= 2:
            return "Middle "
        else:
            return "Highest "

    #plt.style.use(['science', 'ieee'])

    # Load data
    df_avg_all_I0 = pd.read_csv(no_shock_income_path)
    df_avg_all_W0 = pd.read_csv(no_shock_wealth_path)
    df_avg_all_I_plus = pd.read_csv(up_shock_income_path)
    df_avg_all_W_plus = pd.read_csv(up_shock_wealth_path)
    df_avg_all_I_minus = pd.read_csv(down_shock_income_path)
    df_avg_all_W_minus = pd.read_csv(down_shock_wealth_path)

    # Compute the share of overall income and wealth for each percentile
    for df in [df_avg_all_I_plus, df_avg_all_W_plus, df_avg_all_I_minus, df_avg_all_W_minus, df_avg_all_I0, df_avg_all_W0]:
        df['sum'] = df.sum(axis=1)
        df[:] = df.div(df['sum'], axis=0)  # Properly assign the result back to the DataFrame
        df.drop(columns=['sum'], inplace=True)

    # Divide each column by the corresponding column in df_avg_all_I0 (relative change)
    df_avg_all_I_plus = df_avg_all_I_plus.div(df_avg_all_I0) - 1
    df_avg_all_W_plus = df_avg_all_W_plus.div(df_avg_all_W0) - 1
    df_avg_all_I_minus = df_avg_all_I_minus.div(df_avg_all_I0) - 1
    df_avg_all_W_minus = df_avg_all_W_minus.div(df_avg_all_W0) - 1

    # Smoothing if needed
    if smooth:
        df_avg_all_I_plus = df_avg_all_I_plus.rolling(window=window_size).mean()
        df_avg_all_I_minus = df_avg_all_I_minus.rolling(window=window_size).mean()
        df_avg_all_W_plus = df_avg_all_W_plus.rolling(window=window_size).mean()
        df_avg_all_W_minus = df_avg_all_W_minus.rolling(window=window_size).mean()

    # Create plot
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(6, 4))

    axes_right_top = axes[0, 1].twinx()
    axes_right_bottom = axes[1, 1].twinx()
    # Restrict the timeline to rows from index 300 onwards
    df_avg_all_I_plus = df_avg_all_I_plus.iloc[300:]
    df_avg_all_W_plus = df_avg_all_W_plus.iloc[300:]
    df_avg_all_I_minus = df_avg_all_I_minus.iloc[300:]
    df_avg_all_W_minus = df_avg_all_W_minus.iloc[300:]

    # Plotting
    for col in df_avg_all_I_plus.columns:
        label = get_percentile_group(col)
        axes[0, 0].plot(df_avg_all_I_plus.index, df_avg_all_I_plus[col], color=vmap(float(col[11]) / 4), label=label)
        axes[0, 0].axvspan(60, 80, color='gray', alpha=0.1)

    for col in df_avg_all_W_plus.columns:
        axes[1, 0].plot(df_avg_all_W_plus.index, df_avg_all_W_plus[col], color=vmap(float(col[11]) / 4), label=f"{col} ")
        axes[1, 0].axvspan(60, 80, color='gray', alpha=0.1)

    for col in df_avg_all_I_minus.columns:
        axes_right_top.plot(df_avg_all_I_minus.index, df_avg_all_I_minus[col], color=vmap(float(col[11]) / 4))
        axes[0, 1].axvspan(60, 80, color='gray', alpha=0.1)

    for col in df_avg_all_W_minus.columns:
        axes_right_bottom.plot(df_avg_all_W_minus.index, df_avg_all_W_minus[col], color=vmap(float(col[11]) / 4))
        axes[1, 1].axvspan(60, 80, color='gray', alpha=0.1)

    # Set x-axis limits to start at 300
    for ax in [axes[0, 0], axes[1, 0], axes[0, 1], axes[1, 1]]:
        ax.set_xlim(300, None)  # Set minimum x-axis value to 300

    for ax in axes.flatten():
        ax.legend(loc='right', bbox_to_anchor=(1, 0.5))
        axes[0, 0].legend(loc='upper left', fontsize=8)
        axes[0, 0].set_title(r'Fuel price shock +50\%')
        axes[0, 1].set_title(r'Fuel price shock -50\%')
        axes[0, 0].set_ylabel('Relative (real) Income to no shock')
        axes[1, 0].set_ylabel('Relative (real) Wealth to no shock')
        axes[0, 0].get_xaxis().set_visible(False)
        axes[0, 1].get_xaxis().set_visible(False)
        axes[0, 1].get_yaxis().set_visible(False)
        axes[1, 1].get_yaxis().set_visible(False)
        axes[1, 0].set_xlabel('Time')
        axes[1, 1].set_xlabel('Time')
        axes[1, 1].legend().set_visible(False)
        axes[1, 0].legend().set_visible(False)
        axes[0, 1].legend().set_visible(False)
        axes[0, 0].text(0.92, 0.96, 'a)', transform=axes[0, 0].transAxes, fontsize=12, fontweight='bold', va='top')
        axes[0, 1].text(0.92, 0.96, 'b)', transform=axes[0, 1].transAxes, fontsize=12, fontweight='bold', va='top')
        axes[1, 0].text(0.92, 0.96, 'c)', transform=axes[1, 0].transAxes, fontsize=12, fontweight='bold', va='top')
        axes[1, 1].text(0.92, 0.96, 'd)', transform=axes[1, 1].transAxes, fontsize=12, fontweight='bold', va='top')

    fig.patch.set_facecolor('white')
    fig.savefig(output_path, dpi=dpi_value)

    plt.tight_layout()
    #plt.show()







#=======================================================================================================
# Set the working directory to one level down
os.chdir(os.path.join(os.getcwd(), "results"))

plot_new_panel()
plot_KC()

household_data_path = "../data/priceshocks/household_data_down.csv"
cpi_data_path = "../data/priceshocks/down_priceshock.csv"
df_avg_all_I, df_avg_all_W = process_household_data(household_data_path, cpi_data_path)


household_data_path = "../data/priceshocks/household_data_up.csv"
cpi_data_path = "../data/priceshocks/up_priceshock.csv"
df_avg_all_I, df_avg_all_W = process_household_data(household_data_path, cpi_data_path)


household_data_path = "../data/priceshocks/household_data.csv"
cpi_data_path = "../data/priceshocks/no_priceshock.csv"
df_avg_all_I, df_avg_all_W = process_household_data(household_data_path, cpi_data_path)

generate_incomeshock_plot(
    no_shock_income_path="../data/priceshocks/household_data_avg_I.csv",
    no_shock_wealth_path="../data/priceshocks/household_data_avg_W.csv",
    up_shock_income_path="../data/priceshocks/household_data_up_avg_I.csv",
    up_shock_wealth_path="../data/priceshocks/household_data_up_avg_W.csv",
    down_shock_income_path="../data/priceshocks/household_data_down_avg_I.csv",
    down_shock_wealth_path="../data/priceshocks/household_data_down_avg_W.csv",
    output_path="../results/incomeshocks.pdf"
)
