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
    thread_nr::Int64;
    isinput::Bool=false
    )::String

    if isinput
        return folderpath * "gsa_input_run$(run_nr)_thread$(thread_nr).csv"
    else
        return folderpath * "gsa_output_run$(run_nr)_thread$(thread_nr).csv"
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
    n_per_thread::Int64,
    n_threads::Int64,
    inputpath::String;
    samp_strat::String = "lhs"
    )

    # Include Python file containing GSA functions
    @pyinclude("parameters/sensitivity/run_GSA.py")

    for thread_nr in 1:n_threads

        df = DataFrame()

        # Call SAMP function to get the input parameters for running the models
        X = py"call_AAT_sampling"(samp_strat, X_labels, n_per_thread)

        # Add simulation number
        df[!, :sim_nr] = (thread_nr - 1) * n_per_thread + 1 : thread_nr * n_per_thread

        # Add samples parameter values
        for (j,label) in enumerate(keys(X_labels))
            df[!, Symbol(label)] = X[:,j]
        end

        # Write to csv
        filepath = getfilepath(inputpath, run_nr, thread_nr; isinput=true)
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
function generate_simdata_multithreaded(
    X_labels::Dict,
    n_threads::Int64,
    n_per_epoch::Int64,
    n_per_thread::Int64,
    inputpath::String,
    outputpath::String,
    run_nr::Int64;
    t_warmup::Int64=100
    )
    
    Threads.@threads for _ in 1:n_threads

        params = collect(keys(X_labels))
        changedparams = Dict(params .=> 0.0)

        thread_nr = Threads.threadid()

        inputfilepath = getfilepath(inputpath, run_nr, thread_nr; isinput=true)
        outputfilepath = getfilepath(outputpath, run_nr, thread_nr)
        res = nothing

        for (i, row) in enumerate(CSV.Rows(inputfilepath))

            sim_nr = parse(Int64, row.sim_nr)

            # println("sim nr $sim_nr")

            # Fill in changed parameters
            for param in params
                changedparams[param] = parse(Float64, getproperty(row, Symbol(param)))
            end

            # Run the model with changed parameters
            runoutput = run_simulation(
                changed_params=changedparams,
                full_output=false;
                threadnr=Threads.threadid()
            )

            if res == nothing
                res = convertrunoutput(runoutput, sim_nr; return_as_df=true)
            else
                push!(res, (convertrunoutput(runoutput, sim_nr)))
            end

            # If end of epoch is reached, write results to output csv
            if nrow(res) == n_per_epoch || nrow(res) == n_per_thread
                CSV.write(outputfilepath, res; append=i≠n_per_epoch)
                res = nothing
            end
        end
    end
end

getinputpath(threadi::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/input_data/gsa_input_run$(run_nr)_thread$threadi.csv"

getoutputpath(threadi::Int64, run_nr::Int64) = "parameters/sensitivity/sensitivity_runs/output_data/gsa_output_run$(run_nr)_thread$threadi.csv"

"""
Calls SAFE toolbox in Python script.
"""
function run_PAWN(
    X_labels::Dict,
    run_nr::Int64;
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

    labels = collect(keys(X_labels))
    labelnames = ["ψ_E", "μ_1", "κ_{upper}", "ω", "ϵ", "α_{cp}", "p_f", "prog"]
    X = zeros(nrow(df), length(labels))
    X .= df[:, labels]

    # GDP
    # GDP_1st = 100 .* df[!, Symbol("GDP_1st")]
    # py"run_PAWN"(labelnames, X, GDP_1st, "GDP_1", run_nr, "mean GDP growth")

    # GDP_2nd = 100 .* df[!, Symbol("GDP_2nd")]
    # py"run_PAWN"(labelnames, X, GDP_2nd, "GDP_2", run_nr, "var GDP growth")

    # GDP_3rd = 100 .* df[!, Symbol("GDP_3rd")]
    # py"run_PAWN"(labelnames, X, GDP_3rd, "GDP_3", run_nr, "skew GDP growth")

    # GDP_4th = 100 .* df[!, Symbol("GDP_4th")]
    # py"run_PAWN"(labelnames, X, GDP_4th, "GDP_4", run_nr, "kurtosis GDP growth")

    # GDP_acorr = df[!, Symbol("acorr_GDP")]
    # py"run_PAWN"(labelnames, X, GDP_acorr, "acorr_GDP", run_nr, "GDP growth autocorrelation")

    # # GINI coefficients
    # GINI_I_1st = df[!, Symbol("GINI_I_1st")]
    # py"run_PAWN"(labelnames, X, GINI_I_1st, "GINI_I_1", run_nr, "mean income GINI")

    # GINI_I_2nd = df[!, Symbol("GINI_I_2nd")]
    # py"run_PAWN"(labelnames, X, GINI_I_2nd, "GINI_I_2", run_nr, "var income GINI")

    # GINI_W_1st = df[!, Symbol("GINI_W_1st")]
    # py"run_PAWN"(labelnames, X, GINI_W_1st, "GINI_W_1", run_nr, "mean wealth GINI")

    # GINI_W_2nd = df[!, Symbol("GINI_W_2nd")]
    # py"run_PAWN"(labelnames, X, GINI_W_2nd, "GINI_W_2nd", run_nr, "var wealth GINI")
    
    # # Unemployment
    # U_1st = df[!, Symbol("U_1st")]
    # py"run_PAWN"(labelnames, X, U_1st, "U_1", run_nr, "mean unemployment rate")

    # U_2nd = df[!, Symbol("U_2nd")]
    # py"run_PAWN"(labelnames, X, U_2nd, "U_2", run_nr, "var unemployment rate")

    dU_1st = df[!, Symbol("dU1st")]
    py"run_PAWN"(labelnames, X, dU_1st, "dU_1", run_nr, "% change in unemployment rate")

    corr_GDP_dU = df[!, Symbol("corr_GDP_dU")]
    py"run_PAWN"(labelnames, X, corr_GDP_dU, "corr_GDP_dU", run_nr, "corr GDP growth and % change in U")

    # # Poverty
    # FGT_1st = df[!, Symbol("FGT_1st")]
    # py"run_PAWN"(labelnames, X, FGT_1st, "FGT_1", run_nr, "mean FGT")

    # FGT_2nd = df[!, Symbol("FGT_2nd")]
    # py"run_PAWN"(labelnames, X, FGT_2nd, "FGT_2", run_nr, "var FGT")

    # # Carbon emissions
    # em2030 = df[!, Symbol("em2030")]
    # py"run_PAWN"(labelnames, X, em2030, "em2030", run_nr, "CO_2 em 2030")

    # em2040 = df[!, Symbol("em2040")]
    # py"run_PAWN"(labelnames, X, em2040, "em2040", run_nr, "CO_2 em 2040")

    # em2050 = df[!, Symbol("em2050")]
    # py"run_PAWN"(labelnames, X, em2050, "em2050", run_nr, "CO_2 em 2050")

    # # Productivity
    # LP_g_1st = df[!, Symbol("LP_g_1st")]
    # py"run_PAWN"(labelnames, X, LP_g_1st, "LP_g", run_nr, "mean LP growth")

    # EE_g_1st = df[!, Symbol("EE_g_1st")]
    # py"run_PAWN"(labelnames, X, EE_g_1st, "EE_g", run_nr, "mean EE growth")

    # EF_g_1st = df[!, Symbol("EF_g_1st")]
    # py"run_PAWN"(labelnames, X, EF_g_1st, "EF_g", run_nr, "mean EF growth")

end


function parse_commandline(
    M::Int64;
    N_u::Int64=500,
    n::Int64=50,
    N_c::Int64=100
    )

    s = ArgParseSettings()
    @add_arg_table s begin
        
        "--n_sims", "-N"
            help="number of simulations"
            arg_type=Int64
            default=N_u + n * N_c * M
        "--n_per_epoch"
            help="number of simulations per epoch"
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
    end

    return parse_args(s)
end


function main(;
    run_nr::Int64=9,
    n_timesteps::Int64=460,
    n_warmup::Int64=100
    )

    # X_labels = Dict([["α_cp", [0.4, 1.0]],
    #                  ["μ1", [0.0, 0.4]],
    #                  ["ω", [0.0, 1.0]],
    #                  ["ϵ", [0.0, 0.1]],
    #                  ["κ_upper", [0.0, 0.05]],
    #                  ["ψ_E", [0.0, 0.1]],
    #                  ["ψ_Q", [0.0, 0.5]],
    #                  ["ψ_P", [0.0, 0.5]],
    #                  ["p_f", [0.0, 1.0]]]
    #                 )

    X_labels = Dict([
        ["α_cp", [0.6, 1.0]],
        ["prog", [-1.0, 1.0]],
        # ["μ1", [0.0, 0.5]],
        ["ω", [0.0, 1.0]],
        ["λ", [0.0, 1.0]],
        ["ϵ", [0.0, 0.1]],
        ["κ_upper", [0.0, 0.05]],
        # ["ψ_E", [0.0, 0.1]],
        ["p_f", [0.0, 1.0]]]
   )

    parsed_args = parse_commandline(length(X_labels))

    # Unpack parsed arguments
    n_sims = parsed_args["n_sims"]
    n_per_epoch = parsed_args["n_per_epoch"]
    inputpath = parsed_args["inputpath"]
    outputpath = parsed_args["outputpath"]
    geninput = parsed_args["geninput"]
    genoutput = parsed_args["genoutput"]
    runpawn = parsed_args["runpawn"]

    # Specify number of threads and n per thread
    n_threads = Threads.nthreads()
    n_per_thread = ceil(Int64, n_sims / n_threads)

    # Generate parameters used for SA
    if geninput 
        generate_labels(
            X_labels, 
            run_nr,
            n_per_thread,  
            n_threads,
            inputpath
        )
    end

    # Generate simulation data
    if genoutput
        generate_simdata_multithreaded(
            X_labels,
            n_threads,
            n_per_epoch, 
            n_per_thread,
            inputpath,
            outputpath, 
            run_nr
        )
    end

    # Generate PAWN indeces
    if runpawn
        run_PAWN(X_labels, run_nr)
    end
end

main()