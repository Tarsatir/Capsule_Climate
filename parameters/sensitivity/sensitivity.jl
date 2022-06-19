"""
This file contains code used to conduct global sensitivity analysis.
The used method is PAWN (https://www.safetoolbox.info/pawn-method/), 
    which is implemented in Python.
"""

using PyCall

include("../../model/main.jl")

"""
    generate_labels(X_labels::Dict, path::String; N::Int)

Generates labels (parameters) used to run the sensitivity analysis. Saves to csv

    - `X_labels` dict of independent parameters that will be tested.
    - `run_nr`: number of complete run of SA.
    - `N`: number of simulations. 
"""
function generate_labels(
    X_labels::Dict,
    run_nr::Int,
    N::Int,
    N_per_thread::Int,
    n_threads::Int,
    )

    # Define parameter distributions
    samp_strat = "lhs"

    # Call SAMP function to get the input parameters for running the models
    X = py"call_AAT_sampling"(samp_strat, M, X_labels, N)

    for thread_nr in 1:n_threads

        # Define boundaries computed within thread
        thread_start = (thread_nr-1) * N_per_thread + 1
        thread_end = min(thread_nr * N_per_thread, N)

        params = Dict(x=>X[thread_start:thread_end, j] for (j,x) in enumerate(keys(X_labels)))
        
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

        # Add poverty moments
        # TODO

        # Write to csv
        CSV.write(path(run_nr, thread_nr), df)
    end
end


"""
    generate_simdata(X_labels::Dict, n_per_epoch::Int, N_per_thread::Int, run_nr::Int)

Generates data used by sensitivity analysis.
    Code augmented from SAFE package preparation of variables.

    - `X_labels`: dictionary of changed parameters and their ranges.
    - `n_per_epoch`: number of simulations per epoch, after which (intermediate) save is made.
    - `N_per_thread`: number of simulations per thread.
    - `run_nr`: identifier of run attempt.
"""
function generate_simdata(
    X_labels::Dict,
    n_per_epoch::Int,
    N_per_thread::Int,
    run_nr::Int
    )
    
    Threads.@threads for _ in 1:n_threads

        total_completed_runs = 0

        params = collect(keys(X_labels))
        changedparams = Dict(params .=> 0.0)

        n_epochs = ceil(Int64, N_per_thread / n_per_epoch)

        for epoch in 1:n_epochs

            println("Thread $(Threads.threadid()), epoch $epoch")

            df = DataFrame(CSV.File(path(run_nr, Threads.threadid())))

            firstrow = total_completed_runs + 1
            lastrow = firstrow + n_per_epoch

            for row in firstrow:lastrow

                println("   Tr $(Threads.threadid()), row $row")

                # Fill in changed parameters
                for param in params
                    changedparams[param] = df[row, param]
                end

                # Run the model with changed parameters
                GDP_g, GINI_I, GINI_W, U = run_simulation(
                    changed_params=changedparams,
                    full_output=false;
                    threadnr=Threads.threadid()
                )

                # Write GDP data to dataframe
                df[row, "GDP_1st"] = mean(GDP_g)
                df[row, "GDP_2nd"] = var(GDP_g)
                df[row, "GDP_3rd"] = skewness(GDP_g)
                df[row, "GDP_4th"] = kurtosis(GDP_g)

                # Write Gini data to dataframe
                df[row, "Gini_I_1st"] = mean(GINI_I)
                df[row, "Gini_I_2nd"] = var(GINI_I)

                df[row, "Gini_W_1st"] = mean(GINI_W)
                df[row, "Gini_W_2nd"] = var(GINI_W)

                # Write unemployment data to dataframe
                df[row, "U_1st"] = mean(U)
                df[row, "U_2nd"] = var(U)

                # Update total completed runs
                total_completed_runs += 1
                if total_completed_runs == nrow(df)
                    break
                end
            end

            CSV.write(path(run_nr, Threads.threadid()), df)

            if total_completed_runs == nrow(df)
                break
            end
        end
    end

end


"""
Calls SAFE toolbox in Python script.
"""
function run_PAWN(
    X_labels::Dict,
    path::String,
    run_nr;
    N=100
    )

    # Read simulation data, save X and Y as matrices
    X = zeros(N, length(X_labels))
    df = DataFrame(CSV.File(path))

    for (i,label) in enumerate(keys(X_labels))
        X[:,i] = df[!, Symbol(label)]
    end

    labels = collect(keys(X_labels))

    Yg_mean = 100 .* df[!, Symbol("Yg_mean")]
    py"run_PAWN"(labels, X, Yg_mean, "Yg_mean", run_nr, "mean GDP growth")

    Yg_std = df[!, Symbol("Yg_std")]
    py"run_PAWN"(labels, X, Yg_std, "Yg_std", run_nr, "std GDP growth")

    Cg_mean = df[!, Symbol("Cg_mean")]
    py"run_PAWN"(labels, X, Cg_mean, "Cg_mean", run_nr, "mean C growth")

    Cg_std = df[!, Symbol("Cg_std")]
    py"run_PAWN"(labels, X, Cg_std, "Cg_std", run_nr, "std C growth")
end


# Include Python file containing GSA functions
@pyinclude("parameters/sensitivity/run_GSA.py")

run_nr = 5

path(run_nr, thread_nr) = "parameters/sensitivity/sensitivity_runs/sensitivity_run_$(run_nr)_thr_$(thread_nr).csv"

# n_threads = Threads.nthreads()
n_threads = 32

# Define number of similated runs
N_u = 500
n = 50
N_c = 100


# Define number of time steps per simulation and warmup period
n_timesteps = 460
n_warmup = 100

X_labels = Dict([["α_cp", [0.4, 1.0]],
                 ["μ1", [0.0, 0.4]],
                 ["ω", [0.0, 1.0]],
                 ["ϵ", [0.0, 0.1]],
                 ["κ_upper", [0.0, 0.05]],
                 ["ψ_E", [0.0, 0.1]],
                 ["ψ_Q", [0.0, 0.5]],
                 ["ψ_P", [0.0, 0.5]],
                 ["p_f", [0.0, 0.5]]])

# Number of uncertain parameters
M = length(X_labels)

N = N_u + n * N_c * M

# TEMP
# println(N)
N = 35000
N_per_thread = ceil(Int64, N / n_threads)

# Generate parameters used for SA
generate_labels(
    X_labels, 
    run_nr,
    N,
    N_per_thread,  
    n_threads,
)

# Generate simulation data
# n_per_epoch = 2
# generate_simdata(X_labels, n_per_epoch, N_per_thread, run_nr)

# run_PAWN(X_labels, path, run_nr; N=N)