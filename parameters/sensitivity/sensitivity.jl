"""
This file contains code used to conduct global sensitivity analysis.
The used method is PAWN (https://www.safetoolbox.info/pawn-method/), 
    which is implemented in Python.
"""

using PyCall
using ArgParse

include("../../model/main.jl")
include("../../results/experiments/experiment_helpers.jl")


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



"""
    generate_labels(X_labels::Dict, run_nr::Int64, n_per_thread::Int64, 
                    n_threads::Int64, inputpath::String; samp_strat::String = "lhs")

Generates labels (parameters) used to run the sensitivity analysis. Saves to csv

    - `X_labels` dict of independent parameters that will be tested.
    - `run_nr`: number of complete run of SA.
    - `N`: number of simulations.
    - `n_per_thread`: number of simulations per thread
    - `n_threads`: total number of threads
    - `samp_strat` [OPT]: samling strategy (default=`lhs`)
"""
function generate_labels(
    X_labels::Dict,
    run_nr::Int64,
    n_per_parl::Int64,
    n_parl::Int64,
    inputpath::String;
    samp_strat::String = "lhs"
    )

    # Include Python file containing GSA functions
    @pyinclude("parameters/sensitivity/run_GSA.py")

    for parl_id in 1:n_parl

        df = DataFrame()

        # Call SAMP function to get the input parameters for running the models
        X = py"call_AAT_sampling"(samp_strat, X_labels, n_per_parl)

        # Add simulation number
        df[!, :sim_nr] = (parl_id - 1) * n_per_parl + 1 : parl_id * n_per_parl

        # Add samples parameter values
        for (j,label) in enumerate(keys(X_labels))
            df[!, Symbol(label)] = X[:,j]
        end

        # Write to csv
        filepath = getfilepath(inputpath, run_nr, parl_id; isinput=true)
        CSV.write(filepath, df)
    end
end


"""
    generate_simdata(X_labels::Dict, n_per_epoch::Int64, n_per_thread::Int64, run_nr::Int64)

Generates data used by sensitivity analysis.
    Code augmented from SAFE package preparation of variables.

    - `X_labels`: dictionary of changed parameters and their ranges.
    - `n_per_epoch`: number of simulations per epoch, after which (intermediate) save is made.
    - `n_per_thread`: number of simulations per thread.
    - `run_nr`: identifier of run attempt.
"""
function generate_simdata(
    X_labels::Dict,
    Y_labels::Vector{Symbol},
    n_timesteps::Int64,
    n_per_epoch::Int64,
    n_per_parl::Int64,
    t_warmup::Int64,
    inputpath::String,
    outputpath::String,
    sim_nr_only::Int64,
    run_nr::Int64;
    save_full_output::Bool=true,
    proc_nr::Union{Int64, Nothing}=nothing
)

    params = collect(keys(X_labels))
    changedparams = Dict(params .=> 0.0)

    # If multithreaded, take thread id, otherwise take passed sim nr
    parl_nr = isnothing(proc_nr) ? Threads.threadid() : proc_nr

    inputfilepath = getfilepath(inputpath, run_nr, parl_nr; isinput=true)
    outputfilepath = getfilepath(outputpath, run_nr, parl_nr)

    for (i, row) in enumerate(CSV.Rows(inputfilepath))

        # Set seed to the same value for each iteration
        seed = 1234
        Random.seed!(seed)

        sim_nr = parse(Int64, row.sim_nr)

        if sim_nr_only != 0 && sim_nr != sim_nr_only
            continue
        end

        # Fill in changed parameters
        for param in params
            changedparams[param] = parse(Float64, getproperty(row, Symbol(param)))
        end

        # Run the model with changed parameters
        _, model_df = run_simulation(
            T = n_timesteps,
            changed_params = changedparams;
            thread_nr = parl_nr,
            sim_nr = sim_nr
        )


        # Check to save full output per run, or preprocess results first
        if save_full_output

            # Save full time series of selected Y labels
            outputfilepath = get_output_path(sim_nr)
            CSV.write(outputfilepath, model_df[:,Y_labels])


        # TODO: MAKE FUNCTION FOR CONDENSED DATA
        # else
        #     # Preprocess output results
        #     if isnothing(res)
        #         res = convertrunoutput(model_df, sim_nr; return_as_df=true, t_warmup=t_warmup)
        #     else
        #         push!(res, (convertrunoutput(model_df, sim_nr; t_warmup=t_warmup)))
        #     end

        #     # If end of epoch is reached, write results to output csv
        #     if nrow(res) == n_per_epoch || i == n_per_parl
        #         CSV.write(outputfilepath, res; append=i≠n_per_epoch)
        #         res = nothing
        #     end
        end

        # If max number of simulations is reached, stop thread/loop
        if i == n_per_parl
            break
        end
    end
end


# get_input_path(parl_id::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/input_data/gsa_input_run$(run_nr)_thread$parl_id.csv"
# get_output_path(parl_id::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/output_data/gsa_output_run$(run_nr)_thread$parl_id.csv"

get_input_path(parl_id::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/input_data/gsa_input_run$(run_nr)_thread$parl_id.csv"
get_output_path(parl_id::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/output_data/gsa_output_run$(run_nr)_thread$parl_id.csv"
get_output_path(sim_nr::Int64) = "parameters/sensitivity/sensitivity_runs/output_data/gsa_output_run_$(sim_nr).csv"


"""
Calls SAFE toolbox in Python script.
"""
function run_PAWN(
    X_labels::Dict,
    run_nr::Int64,
    N_u::Int64,
    N_c::Int64,
    outputpath::String;
    nthreads::Int64=16
)

    # Include Python file containing GSA functions
    @pyinclude("parameters/sensitivity/run_GSA.py")

    # Open all input dataframes
    dfs = [DataFrame(CSV.File(get_input_path(i, run_nr))) for i in 1:nthreads]
    input_df = dfs[1]
    for df in dfs[2:end]
        append!(input_df, df)
    end

    # Open all output dataframes
    dfs = [DataFrame(CSV.File(get_output_path(i, run_nr))) for i in 1:nthreads]
    output_df = dfs[1]
    for df in dfs[2:end]
        append!(output_df, df)
    end

    # Merge input and output dataframes
    df = innerjoin(input_df, output_df, on=:sim_nr)

    # Define labels and label names
    labels = ["α_cp", "ω", "λ", "κ_upper", "prog", "μ1", "ϵ_w", "ϵ_μ", "p_f", "ψ_P", "ψ_Q", "ψ_E"]
    labelnames = ["\\alpha_{cp}", "\\omega", "\\lambda", "\\chi_{upper}", "prog", "\\mu_1", 
                  "\\epsilon_w", "\\bar{\\epsilon}_\\mu", "p_f", "\\psi_P", "\\psi_Q", "\\psi_E"]

    X = zeros(nrow(df), length(labels))
    X .= df[:, labels]

    # Compute the critical value
    crit = 1.358 * sqrt((N_u + N_c) / (N_u * N_c))

    # Make dataframe to collect KS values
    KS_labels = vcat(["y"], labels)
    df_KS = DataFrame(vcat([String[]], [Float64[] for _ in KS_labels[2:end]]), KS_labels)

    Y = [
        ["GDP_1st", "mean \$\\Delta GDP\$"],
        ["GDP_2nd", "var \$\\Delta GDP\$"],
        ["GDP_3rd", "skew \$\\Delta GDP\$"],
        ["GDP_4th", "kurtosis \$\\Delta GDP\$"],
        ["dQ_1st", "mean \$\\Delta Q\$"],
        ["dQ_2nd", "var \$\\Delta Q\$"],
        ["LIS_1st", "mean LIS"],
        ["GINI_I_1st", "mean income GINI"],
        ["GINI_I_2nd", "var income GINI"],
        ["GINI_W_1st", "mean wealth GINI"],
        ["GINI_W_2nd", "var wealth GINI"],
        ["U_1st", "mean U"],
        ["U_2nd", "var U"],
        ["dI_1st", "mean \$\\Delta I\$"],
        ["dC_1st", "mean \$\\Delta C\$"],
        ["bankr_1st", "mean bankr"],
        ["LP_g_1st", "mean \$\\Delta LP\$"],
        ["EE_g_1st", "mean \$\\Delta EE\$"],
        ["EF_g_1st", "mean \$\\Delta LP\$"],
        ["em2030", "\$CO_2\$ index 2030"],
        ["em2040", "\$CO_2\$ index 2040"],
        ["em2050", "\$CO_2\$ index 2050"]
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


function parse_commandline(
    M::Int64,
    N_u::Int64,
    n::Int64,
    N_c::Int64
)

    s = ArgParseSettings()
    @add_arg_table s begin
        
        "--n_sims", "-N"
            help="number of simulations"
            arg_type=Int64
            default=N_u + n * N_c * M
        "--n_per_epoch"
            help="number of simulations per epoch (after which data is saved)"
            arg_type=Int64
            default=50
        "--inputpath"
            help="path to directory with input files"
            arg_type=String
            default="parameters/sensitivity/sensitivity_runs/input_data/"
        "--outputpath"
            help="path to directory with output files"
            arg_type=String
            default="parameters/sensitivity/sensitivity_runs/output_data/"
        "--n_processes"
            help="total number of processes"
            arg_type=Int64
            default=0
        "--process_id"
            help="process id, if passed, program is assumed to be executed multiprocessing, 
                  otherwise assumed to be multithreaded."
            arg_type=Int64
            default=0
        "--geninput"
            help="determines if input parameters are generated"
            arg_type=Bool
            default=false
        "--genoutput"
            help="determines if output parameters are generated (may take a long time)"
            arg_type=Bool
            default=true
        "--runpawn"
            help="determines if PAWN indeces are run"
            arg_type=Bool
            default=false
        "--sim_nr"
            help="run a single simulation (for testing purposes), sim_nr must correspond with process_id"
            arg_type=Int64
            default=0
    end

    return parse_args(s)
end


function main(;
    run_nr::Int64=9,
    n_timesteps::Int64=10,
    t_warmup::Int64=10
)

    # Set parameters for optimal sampling quantity
    N_u = 500
    n = 20
    N_c = 100

    # Define parameters that are varied over and their ranges
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

    # Dependent variables that are saved  !!! NOTE: Change SA output here! !!!
    Y_labels = [:GDP, :emissions_index, :energy_percentage ]

    parsed_args = parse_commandline(length(X_labels), N_u, n, N_c)

    # Unpack parsed arguments
    n_sims = parsed_args["n_sims"]
    n_per_epoch = parsed_args["n_per_epoch"]
    inputpath = parsed_args["inputpath"]
    outputpath = parsed_args["outputpath"]
    n_processes = parsed_args["n_processes"]
    process_id = parsed_args["process_id"]
    geninput = parsed_args["geninput"]
    genoutput = parsed_args["genoutput"]
    runpawn = parsed_args["runpawn"]
    sim_nr_only = parsed_args["sim_nr"]

    
    if n_processes == 0
        # Program assumed to be multithreaded
        n_threads = Threads.nthreads()
        n_per_parl = ceil(Int64, n_sims / n_threads)
    else
        n_per_parl = ceil(Int64, n_sims / n_processes)
    end

    n_parl = n_processes == 0 ? n_threads : n_processes

    # Generate parameters used for SA
    if geninput 
        generate_labels(
            X_labels, 
            run_nr,
            n_per_parl,  
            n_parl,
            inputpath
        )
    end

    # Generate simulation data
    if genoutput
        if n_processes == 0
            Threads.@threads for _ in 1:n_parl
                generate_simdata(
                    X_labels,
                    Y_labels,
                    n_timesteps,
                    n_per_epoch, 
                    n_per_parl,
                    t_warmup,
                    inputpath,
                    outputpath,
                    sim_nr_only, 
                    run_nr
                )
            end
        # else
        #     generate_simdata(
        #         X_labels,
        #         Y_labels,
        #         n_per_epoch, 
        #         n_per_parl,
        #         t_warmup,
        #         inputpath,
        #         outputpath, 
        #         run_nr;
        #         proc_nr=process_id,
        #         sim_nr_only=sim_nr
        #     )
        end
    end

    # Generate PAWN indices
    if runpawn
        run_PAWN(X_labels, run_nr, N_u, N_c, outputpath)
    end
end

main()