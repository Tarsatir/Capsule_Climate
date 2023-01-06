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
    n_per_epoch::Int64,
    n_per_parl::Int64,
    t_warmup::Int64,
    inputpath::String,
    outputpath::String,
    run_nr::Int64;
    proc_nr::Union{Int64, Nothing}=nothing,
    sim_nr_only::Int64
)

    params = collect(keys(X_labels))
    changedparams = Dict(params .=> 0.0)

    # If multithreaded, take thread id, otherwise take passed sim nr
    parl_nr = proc_nr == nothing ? Threads.threadid() : proc_nr

    inputfilepath = getfilepath(inputpath, run_nr, parl_nr; isinput=true)
    outputfilepath = getfilepath(outputpath, run_nr, parl_nr)
    res = nothing

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
        runoutput = run_simulation(
            changed_params = changedparams;
            threadnr = parl_nr,
            sim_nr = sim_nr
        )

        if isnothing(res)
            res = convertrunoutput(runoutput, sim_nr; return_as_df=true, t_warmup=t_warmup)
        else
            push!(res, (convertrunoutput(runoutput, sim_nr; t_warmup=t_warmup)))
        end

        # If end of epoch is reached, write results to output csv
        if nrow(res) == n_per_epoch || i == n_per_parl
            CSV.write(outputfilepath, res; append=i≠n_per_epoch)
            res = nothing
        end

        if i == n_per_parl
            break
        end
    end
end


getinputpath(parl_id::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/input_data/gsa_input_run$(run_nr)_thread$parl_id.csv"
getoutputpath(parl_id::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/output_data/gsa_output_run$(run_nr)_thread$parl_id.csv"


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
    dfs = [DataFrame(CSV.File(getinputpath(i, run_nr))) for i in 1:nthreads]
    input_df = dfs[1]
    for df in dfs[2:end]
        append!(input_df, df)
    end

    # Open all output dataframes
    dfs = [DataFrame(CSV.File(getoutputpath(i, run_nr))) for i in 1:nthreads]
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

    for (y, y_title) in Y
        ts = df[!, Symbol(y)]
        ts_KS = py"run_PAWN"(labelnames, X, ts, y, run_nr, y_title, crit)
        push!(df_KS, vcat([y], ts_KS))
        CSV.write(outputfilepath, df_KS)
    end

    # GDP
    # GDP_1st = 100 .* df[!, Symbol("GDP_1st")]
    # GDP_1_KS = py"run_PAWN"(labelnames, X, GDP_1st, "GDP_1", run_nr, "mean \$\\Delta GDP\$", crit)
    # push!(df_KS, vcat(["GDP_1st"], GDP_1_KS))

    # GDP_2nd = 100 .* df[!, Symbol("GDP_2nd")]
    # GDP_2_KS = py"run_PAWN"(labelnames, X, GDP_2nd, "GDP_2", run_nr, "var \$\\Delta GDP\$", crit)
    # push!(df_KS, vcat(["GDP_2nd"], GDP_2_KS))

    # GDP_3rd = 100 .* df[!, Symbol("GDP_3rd")]
    # GDP_3_KS = py"run_PAWN"(labelnames, X, GDP_3rd, "GDP_3", run_nr, "skew \$\\Delta GDP\$", crit)
    # push!(df_KS, vcat(["GDP_3rd"], GDP_3_KS))

    # GDP_4th = 100 .* df[!, Symbol("GDP_4th")]
    # GDP_4_KS = py"run_PAWN"(labelnames, X, GDP_4th, "GDP_4", run_nr, "kurtosis \$\\Delta GDP\$", crit)
    # push!(df_KS, vcat(["GDP_4th"], GDP_4_KS))

    # GDP_acorr = df[!, Symbol("acorr_GDP")]
    # GDP_ac_KS = py"run_PAWN"(labelnames, X, GDP_acorr, "acorr_GDP", run_nr, "\$\\Delta GDP\$ autocorrelation", crit)
    # push!(df_KS, vcat(["GDP_ac"], GDP_ac_KS))

    # dQ_1st = 100 .* df[!, Symbol("dQ_1st")]
    # dQ_1_KS = py"run_PAWN"(labelnames, X, dQ_1st, "dQ_1", run_nr, "mean \$\\Delta Q\$", crit)
    # push!(df_KS, vcat(["dQ_1st"], dQ_1_KS))

    # LIS_1st = df[!, Symbol("LIS")]

    # GINI coefficients
    # GINI_I_1st = df[!, Symbol("GINI_I_1st")]
    # GINI_I_1_KS = py"run_PAWN"(labelnames, X, GINI_I_1st, "GINI_I_1", run_nr, "mean income GINI", crit)
    # push!(df_KS, vcat(["GINI_I_1"], GINI_I_1_KS))

    # GINI_I_2nd = df[!, Symbol("GINI_I_2nd")]
    # GINI_I_2_KS = py"run_PAWN"(labelnames, X, GINI_I_2nd, "GINI_I_2", run_nr, "var income GINI", crit)
    # push!(df_KS, vcat(["GINI_I_2"], GINI_I_2_KS))

    # GINI_W_1st = df[!, Symbol("GINI_W_1st")]
    # GINI_W_1_KS = py"run_PAWN"(labelnames, X, GINI_W_1st, "GINI_W_1", run_nr, "mean wealth GINI", crit)
    # push!(df_KS, vcat(["GINI_W_1"], GINI_W_1_KS))

    # GINI_W_2nd = df[!, Symbol("GINI_W_2nd")]
    # py"run_PAWN"(labelnames, X, GINI_W_2nd, "GINI_W_2nd", run_nr, "var wealth GINI", crit)
    
    # Unemployment
    # U_1st = df[!, Symbol("U_1st")]
    # U_1_KS = py"run_PAWN"(labelnames, X, U_1st, "U_1", run_nr, "mean U", crit)
    # push!(df_KS, vcat(["U_1"], U_1_KS))

    # U_2nd = df[!, Symbol("U_2nd")]
    # U_2_KS = py"run_PAWN"(labelnames, X, U_2nd, "U_2", run_nr, "var U", crit)
    # push!(df_KS, vcat(["U_2"], U_2_KS))

    # dU_1st = df[!, Symbol("dU1st")]
    # py"run_PAWN"(labelnames, X, dU_1st, "dU_1", run_nr, "\$\\Delta U\$")

    # corr_GDP_dU = df[!, Symbol("corr_GDP_dU")]
    # py"run_PAWN"(labelnames, X, corr_GDP_dU, "corr_GDP_dU", run_nr, "corr \$\\Delta GDP\$ and \$\\Delta U\$")

    # Poverty
    # FGT_1st = df[!, Symbol("FGT_1st")]
    # FGT_1_KS = py"run_PAWN"(labelnames, X, FGT_1st, "FGT_1", run_nr, "mean FGT", crit)
    # push!(df_KS, vcat(["FTG_1"], FGT_1_KS))

    # FGT_2nd = df[!, Symbol("FGT_2nd")]
    # py"run_PAWN"(labelnames, X, FGT_2nd, "FGT_2", run_nr, "var FGT")

    # # Carbon emissions
    # em2030 = df[!, Symbol("em2030")]
    # em2030_KS = py"run_PAWN"(labelnames, X, em2030, "em2030", run_nr, "\$CO_2\$ em 2030", crit)
    # push!(df_KS, vcat(["em2030"], em2030_KS))

    # em2040 = df[!, Symbol("em2040")]
    # em2040_KS = py"run_PAWN"(labelnames, X, em2040, "em2040", run_nr, "\$CO_2\$ em 2040", crit)
    # push!(df_KS, vcat(["em2040"], em2040_KS))

    # em2050 = df[!, Symbol("em2050")]
    # em2050_KS = py"run_PAWN"(labelnames, X, em2050, "em2050", run_nr, "\$CO_2\$ em 2050", crit)
    # push!(df_KS, vcat(["em2050"], em2050_KS))

    # # # Productivity
    # LP_g_1st = df[!, Symbol("LP_g_1st")]
    # LP_1_KS = py"run_PAWN"(labelnames, X, LP_g_1st, "LP_g", run_nr, "mean LP growth", crit)
    # push!(df_KS, vcat(["LP_1"], LP_1_KS))

    # EE_g_1st = df[!, Symbol("EE_g_1st")]
    # EE_1_KS = py"run_PAWN"(labelnames, X, EE_g_1st, "EE_g", run_nr, "mean EE growth", crit)
    # push!(df_KS, vcat(["EE_1"], EE_1_KS))

    # EF_g_1st = df[!, Symbol("EF_g_1st")]
    # EF_1_KS = py"run_PAWN"(labelnames, X, EF_g_1st, "EF_g", run_nr, "mean EF growth", crit)
    # push!(df_KS, vcat(["EF_1"], EF_1_KS))
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
    n_timesteps::Int64=660,
    t_warmup::Int64=300
)

    # Set parameters for optimal sampling quantity
    N_u = 500
    n = 20
    N_c = 100

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
    sim_nr = parsed_args["sim_nr"]

    
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
                    n_per_epoch, 
                    n_per_parl,
                    t_warmup,
                    inputpath,
                    outputpath, 
                    run_nr
                )
            end
        else
            generate_simdata(
                X_labels,
                n_per_epoch, 
                n_per_parl,
                t_warmup,
                inputpath,
                outputpath, 
                run_nr;
                proc_nr=process_id,
                sim_nr_only=sim_nr
            )
        end
    end

    # Generate PAWN indices
    if runpawn
        run_PAWN(X_labels, run_nr, N_u, N_c, outputpath)
    end
end

main()