ENV["PYTHON"] = "/home/imengesha/Climate_Paper/Climate-Paper/myenv/bin/python"
using PyCall
using ArgParse
using DataFrames

include("../../model/main.jl")
include("../../results/experiments/helpers.jl")
#include("sensitivity.jl")

# Include necessary packages
using DataFrames
using CSV
# ... any other required packages

# Include other required Julia files
# include("other_required_script.jl")

function getfilepath(
    folderpath::String,
    run_nr::Int64,
    parl_id::Int64;
    isinput::Bool=false
    )::String

    if isinput
        return folderpath * "gsa_input_run$(run_nr)_thread$(parl_id).csv"
    else
        return folderpath * "gsa_output_run$(run_nr)_thread$(parl_id).csv"
    end
end

get_input_path1(parl_id::Int64, run_nr::Int64) = "sensitivity_runs/input_data/gsa_input_run$(run_nr)_thread$parl_id.csv"
get_output_path1(parl_id::Int64, run_nr::Int64) = "sensitivity_runs/output_data/gsa_output_run$(run_nr)_thread$parl_id.csv"

function run_PAWN1(
    X_labels::Dict,
    run_nr::Int64,
    N_u::Int64,
    N_c::Int64,
    outputpath::String;
    nthreads::Int64=8
)

    # Include Python file containing GSA functions
    @pyinclude("run_GSA.py")

    # Open all input dataframes
    dfs = [DataFrame(CSV.File(get_input_path1(i, run_nr))) for i in 1:nthreads]
    input_df = dfs[1]
    for df in dfs[2:end]
        append!(input_df, df)
    end

    # Open all output dataframes
    dfs = [DataFrame(CSV.File(get_output_path1(i, run_nr))) for i in 1:nthreads]
    output_df = dfs[1]
    for df in dfs[2:end]
        append!(output_df, df)
    end

    # Merge input and output dataframes
    df = DataFrames.innerjoin(input_df, output_df, on=:sim_nr)

    # Define labels and label names
    labels = [
    "α_maxdev",
    "ω",
    "λ",
    "κ_upper",
    "prog",
    "μ1",
    "ϵ_w",
    "ϵ_μ",
    "p_f",
    "ψ_P",
    "ψ_Q",
    "ψ_E"
    ]

    labelnames = [
        "\\alpha_{\\text{maxdev}}",
        "\\omega",
        "\\lambda",
        "\\kappa_{\\text{upper}}",
        "\\text{prog}",
        "\\mu_1",
        "\\epsilon_w",
        "\\epsilon_\\mu",
        "p_f",
        "\\psi_P",
        "\\psi_Q",
        "\\psi_E"
    ]

    X = zeros(nrow(df), length(labels))
    X .= df[:, labels]

    # Compute the critical value
    crit = 1.358 * sqrt((N_u + N_c) / (N_u * N_c))

    # Make dataframe to collect KS values
    KS_labels = vcat(["y"], labels)
    df_KS = DataFrame(vcat([String[]], [Float64[] for _ in KS_labels[2:end]]), KS_labels)

    Y = [
        ["em_index", "CO_2 Emission Index"],
        ["energy_percentage", "Energy Percentage of Total"],
        ["GDP", "Gross Domestic Product"],
        ["bankrupt_cp", "Bankrupties of CP"],
        ["GINI_I", "PGINI_I"],
        ["GINI_W", "GINI_W"],
        ["LIS", "Labor Income Share"],
        ["avg_pi_LP", "avg_pi_LP"],
        ["avg_A_EE", "avg_A_EE"],
        ["avg_A_EF", "avg_A_EF"],
        ["U", "Unemployment"],
        ["total_I", "Investment"],
        ["total_C", "Consumption"]
    ]
    
    

    # Write results to csv
    outputfilepath = outputpath * "KS_mean_values.csv"

    # Run PAWN for all input parameters
    for (y, y_title) in Y
        ts = df[!, Symbol(y)]
        ts_KS = py"run_PAWN"(labelnames, X, ts, y, run_nr, y_title, crit)
        push!(df_KS, vcat([y], ts_KS))
        CSV.write(outputfilepath, df_KS)
    end
end

# Define parameters to pass to the function
X_labels = Dict([
    ["α_maxdev", [0.005, 0.5]],
    ["ρ", [0.05, 0.8]],
    ["prog", [-1.0, 1.0]],
    ["μ1", [0.0, 0.5]],
    ["ω", [0.0, 1.0]],
    ["λ", [0.0, 1.0]],
    ["ϵ_w", [0.0, 0.1]],
    ["ϵ_μ", [0.0, 0.1]],
    ["κ_upper", [0.0, 0.01]],
    ["ψ_E", [0., 0.25]],
    ["ψ_Q", [0., 0.25]],
    ["ψ_P", [0., 0.25]],
    ["p_f", [0.0, 1.0]]
])

run_nr = 9  # example run number
N_u = 500
N_c = 100 # number of conditional samples
outputpath = "sensitivity_runs/output_data/"  # replace with your actual path

# Call the function with your parameters
run_PAWN1(X_labels, run_nr, N_u, N_c, outputpath; nthreads=8)
