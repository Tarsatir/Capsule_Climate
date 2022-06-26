"""
This file contains code used to conduct global sensitivity analysis.
The used method is PAWN (https://www.safetoolbox.info/pawn-method/), 
    which is implemented in Python.
"""

using PyCall
using ArgParse

include("../../model/main.jl")

"""
    generate_labels(X_labels::Dict, path::String; N::Int64)

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
    N::Int64,
    n_per_thread::Int64,
    n_threads::Int64;
    samp_strat::String = "lhs"
    )

    for thread_nr in 1:n_threads

        # Call SAMP function to get the input parameters for running the models
        X = py"call_AAT_sampling"(samp_strat, M, X_labels, n_per_thread)

        # Define boundaries computed within thread
        # thread_start = (thread_nr-1) * n_per_thread + 1
        # thread_end = min(thread_nr * n_per_thread, N)

        # params = Dict(x=>X[thread_start:thread_end, j] for (j,x) in enumerate(keys(X_labels)))
        params = Dict(x=>X[:, j] for (j,x) in enumerate(keys(X_labels)))
        
        # Write to dataframe
        df = DataFrame(params)

        # Add GDP moments
        df[!, "GDP_1st"] .= NaN
        df[!, "GDP_2nd"] .= NaN
        df[!, "GDP_3rd"] .= NaN
        df[!, "GDP_4th"] .= NaN

        # Add Gini moments
        df[!, "Gini_I_1st"] .= NaN
        df[!, "Gini_I_2nd"] .= NaN

        df[!, "Gini_W_1st"] .= NaN
        df[!, "Gini_W_2nd"] .= NaN

        # Add unemployment moments
        df[!, "U_1st"] .= NaN
        df[!, "U_2nd"] .= NaN

        # Add productivity growth
        df[!, "LP_1st"] .= NaN
        df[!, "LP_2nd"] .= NaN

        df[!, "EE_1st"] .= NaN
        df[!, "EE_2nd"] .= NaN

        df[!, "EF_1st"] .= NaN
        df[!, "EF_2nd"] .= NaN

        # Add poverty moments
        df[!, "FGT_1st"] .= NaN
        df[!, "FGT_2nd"] .= NaN

        # Add emissions indexes
        df[!, "em2030"] .= NaN
        df[!, "em2040"] .= NaN
        df[!, "em2050"] .= NaN

        # Write to csv
        CSV.write(path(run_nr, thread_nr), df)
    end
end


function compute_growthrates(
    ts::Vector{Float64}
    )

    return 100 .* (ts[1:end-1] .- ts[2:end]) ./ ts[2:end]
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
    n_threads::Int64,
    n_per_epoch::Int64,
    n_per_thread::Int64,
    run_nr::Int64;
    t_warmup::Int64=100
    )
    
    Threads.@threads for _ in 1:n_threads

        total_completed_runs = 0

        params = collect(keys(X_labels))
        changedparams = Dict(params .=> 0.0)

        # n_epochs = ceil(Int64, n_per_thread / n_per_epoch)

        df = DataFrame(CSV.File(path(run_nr, Threads.threadid())))

        # for _ in 1:n_epochs

        for row in 1:n_per_thread

            # df = DataFrame(CSV.File(path(run_nr, Threads.threadid())))

            # Update firstrow and 
            # if total_completed_runs % n_per_epoch == 0
                # firstrow = total_completed_runs + 1
                # lastrow = firstrow + n_per_epoch - 1
            # end

            # println("       $(firstrow), $(lastrow)")

            # for row in firstrow:lastrow

            println("   Tr $(Threads.threadid()), row $row")

            # Fill in changed parameters
            for param in params
                changedparams[param] = df[row, param]
            end

            # Run the model with changed parameters
            runoutput = run_simulation(
                changed_params=changedparams,
                full_output=false;
                threadnr=Threads.threadid()
            )

            # Write GDP data to dataframe
            df[row, "GDP_1st"] = mean(runoutput.GDP_growth[t_warmup:end])
            df[row, "GDP_2nd"] = var(runoutput.GDP_growth[t_warmup:end])
            df[row, "GDP_3rd"] = skewness(runoutput.GDP_growth[t_warmup:end])
            df[row, "GDP_4th"] = kurtosis(runoutput.GDP_growth[t_warmup:end])

            # Write unemployment data to dataframe
            df[row, "U_1st"] = mean(runoutput.U[t_warmup:end])
            df[row, "U_2nd"] = var(runoutput.U[t_warmup:end])

            # Write Gini data to dataframe
            df[row, "Gini_I_1st"] = mean(runoutput.GINI_I[t_warmup:end])
            df[row, "Gini_I_2nd"] = var(runoutput.GINI_I[t_warmup:end])

            df[row, "Gini_W_1st"] = mean(runoutput.GINI_W[t_warmup:end])
            df[row, "Gini_W_2nd"] = var(runoutput.GINI_W[t_warmup:end])

            # Add productivity growth
            LP_g = compute_growthrates(runoutput.avg_π_LP)
            df[row, "LP_1st"] .= mean(LP_g[t_warmup:end])
            df[row, "LP_2nd"] .= var(LP_g[t_warmup:end])

            EE_g = compute_growthrates(runoutput.avg_π_EE)
            df[row, "EE_1st"] .= mean(EE_g[t_warmup:end])
            df[row, "EE_2nd"] .= var(EE_g[t_warmup:end])

            EF_g = compute_growthrates(runoutput.avg_π_EF)
            df[row, "EF_1st"] .= mean(EF_g[t_warmup:end])
            df[row, "EF_2nd"] .= var(EF_g[t_warmup:end])

            # Write poverty data to dataframe
            df[row, "FGT_1st"] = mean(runoutput.FGT[t_warmup:end])
            df[row, "FGT_2nd"] = var(runoutput.FGT[t_warmup:end])

            # Write emissions indexes
            df[row, "em2030"] = runoutput.emissions_index[220]
            df[row, "em2040"] = runoutput.emissions_index[340]
            df[row, "em2050"] = runoutput.emissions_index[460]

            # Update total completed runs
            total_completed_runs += 1

                # if total_completed_runs == nrow(df)
                #     break
                # end
            # end

            if total_completed_runs % n_per_epoch == 0
                println("   writing to csv...")
                # firstrow = total_completed_runs + 1
                # lastrow = firstrow + n_per_epoch - 1
                CSV.write(path(run_nr, Threads.threadid()), df)
            end
            # CSV.write(path(run_nr, Threads.threadid()), df)

            if total_completed_runs == nrow(df)
                return nothing
            end
        end
    end
end


getpath(threadi::Int64) = "parameters/sensitivity/sensitivity_runs_10000/sensitivity_run_6_thr_$threadi.csv"

"""
Calls SAFE toolbox in Python script.
"""
function run_PAWN(
    X_labels::Dict,
    run_nr;
    nthreads::Int64=16,
    N=100
    )

    # Read simulation data, save X and Y as matrices
    # df = DataFrame(CSV.File(path))

    # Open all dataframes
    dfs = [DataFrame(CSV.File(getpath(i))) for i in 1:nthreads]

    all_df = dfs[1]
    for df in dfs[2:end]
        append!(all_df, df)
    end

    filter!(row -> !isnan(row.GDP_1st), all_df)

    labels = collect(keys(X_labels))
    X = zeros(nrow(all_df), length(labels))
    X .= all_df[:, labels]

    # GDP
    # GDP_1st = 100 .* all_df[!, Symbol("GDP_1st")]
    # py"run_PAWN"(labels, X, GDP_1st, "GDP_1", run_nr, "mean GDP growth")

    GDP_2nd = 100 .* all_df[!, Symbol("GDP_2nd")]
    py"run_PAWN"(labels, X, GDP_2nd, "GDP_2", run_nr, "var GDP growth")

    GDP_3rd = 100 .* all_df[!, Symbol("GDP_3rd")]
    py"run_PAWN"(labels, X, GDP_3rd, "GDP_3", run_nr, "skew GDP growth")

    GDP_4th = 100 .* all_df[!, Symbol("GDP_4th")]
    py"run_PAWN"(labels, X, GDP_4th, "GDP_4", run_nr, "kurtosis GDP growth")

    # GINI coefficients
    # GINI_I_1st = all_df[!, Symbol("Gini_I_1st")]
    # py"run_PAWN"(labels, X, GINI_I_1st, "GINI_I_1", run_nr, "mean income GINI")

    # GINI_W_1st = all_df[!, Symbol("Gini_W_1st")]
    # py"run_PAWN"(labels, X, GINI_W_1st, "GINI_W_1", run_nr, "mean wealth GINI")
    
    # Unemployment
    # U_1st = all_df[!, Symbol("U_1st")]
    # py"run_PAWN"(labels, X, U_1st, "U_1", run_nr, "mean unemployment rate")

    # Poverty
    # FGT_1st = all_df[!, Symbol("FGT_1st")]
    # py"run_PAWN"(labels, X, FGT_1st, "FGT_1", run_nr, "mean FGT")

    # Carbon emissions
    # carbon_conc = all_df[!, Symbol("carbon_conc")]
    # py"run_PAWN"(labels, X, carbon_conc, "carbon_conc", run_nr, "carbon concentration")

end

# Define number of similated runs
N_u = 500
n = 50
N_c = 100

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
            default=10
        "--path"
            help="path to directory with parameter and result files"
            arg_type=String
            default="parameters/sensitivity/sensitivity_runs/"
    end

    return parse_args(s)
end


# Include Python file containing GSA functions
# @pyinclude("parameters/sensitivity/run_GSA.py")


function main(;
    run_nr::Int64=6
    )

    X_labels = Dict([["α_cp", [0.4, 1.0]],
                     ["μ1", [0.0, 0.4]],
                     ["ω", [0.0, 1.0]],
                     ["ϵ", [0.0, 0.1]],
                     ["κ_upper", [0.0, 0.05]],
                     ["ψ_E", [0.0, 0.1]],
                     ["ψ_Q", [0.0, 0.5]],
                     ["ψ_P", [0.0, 0.5]],
                     ["p_f", [0.0, 1.0]]]
                    )

    parsed_args = parse_commandline(length(X_labels))
    
    println(parsed_args)

    # if length(ARGS) == 0
    #     dir = "parameters/sensitivity/sensitivity_runs/"
    # else
    #     dir = ARGS[1]
    # end

    # path(run_nr, thread_nr) = "$(dir)sensitivity_run_$(run_nr)_thr_$(thread_nr).csv"



    # # Define number of time steps per simulation and warmup period
    # n_timesteps = 460
    # n_warmup = 100



    # # Number of uncertain parameters
    # M = length(X_labels)

    # if length(ARGS) < 2
    #     N = N_u + n * N_c * M
    # else
    #     N = parse(Int64, ARGS[2])
    # end

    # # Specify number of threads and n per thread
    # n_threads = Threads.nthreads()
    # n_per_thread = ceil(Int64, N / n_threads)

    # # Generate parameters used for SA
    # # generate_labels(
    # #     X_labels, 
    # #     run_nr,
    # #     N,
    # #     n_per_thread,  
    # #     n_threads,
    # # )

    # # Generate simulation data
    # # n_per_epoch = 10
    # n_per_epoch = 2
    # @time generate_simdata(X_labels, n_threads, n_per_epoch, n_per_thread, run_nr)

    # # run_PAWN(X_labels, run_nr; N=N)
end

main()