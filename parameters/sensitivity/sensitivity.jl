"""
This file contains code used to conduct global sensitivity analysis.
The used method is PAWN (https://www.safetoolbox.info/pawn-method/), 
    which is implemented in Python.
"""

using PyCall

include("../../model/main.jl")


"""
Generates data used by sensitivity analysis.
    Code augmented from SAFE package preparation of variables.
"""
function generate_data(
    X_labels::Dict,
    path::String;
    N=100::Int
    )

    # Set up which parameter to apply GSA to, and the range of the parameter

    # Number of uncertain parameters
    M = length(X_labels)

    # Define parameter distributions
    samp_strat = "lhs"

    # Call SAMP function to get the input parameters for running the models
    X = py"call_AAT_sampling"(samp_strat, M, X_labels, N)
    
    # Run simulations for the given parameter values
    Y = zeros(N)
    Yg = zeros(N)

    for i in 1:N

        println("Simulation nr: $i")

        # Retrieve variables, repackage to pass to simulation
        Xi = X[i,:]
        changed_params = Dict(k=>Xi[i] for (i,k) in enumerate(keys(X_labels)))

        # println(changed_params)

        # Run simulation
        Yi, Ygi = run_simulation(
                    changed_params=changed_params,
                    full_output=false
                    )

        # Store output data
        Y[i] = Yi
        Yg[i] = Ygi
    end

    # Set up dataframe containing results
    df = DataFrame(Dict(x=>X[:,i] for (i,x) in enumerate(keys(X_labels))))
    df[!, "GDP"] = Y
    df[!, "GDP_growth"] = Yg

    # Write results to csv
    CSV.write(path, df)

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

    Y = df[!, Symbol("GDP")]
    py"run_PAWN"(collect(keys(X_labels)), X, Y, run_nr, "GDP")

    Yg = df[!, Symbol("GDP_growth")]
    py"run_PAWN"(collect(keys(X_labels)), X, Yg, run_nr, "GDP growth")


end


# Include Python file containing GSA functions
@pyinclude("parameters/sensitivity/run_GSA.py")

run_nr = 3

path = "parameters/sensitivity/sensitivity_runs/sensitivity_run_$(run_nr).csv"

N = 1

X_labels = Dict([["α_cp", [0.6, 1.0]],
                 ["μ1", [0.0, 0.5]],
                 ["ω", [0.0, 1.0]],
                 ["ϵ", [0.01, 0.1]],
                 ["κ_upper", [0.01, 0.3]],
                 ["κ_lower", [-0.3, -0.01]],
                 ["β1", [1.0, 5.0]],
                 ["β2", [1.0, 5.0]],
                 ["ψ_E", [0.0, 1.0]],
                 ["ψ_Q", [0.0, 1.0]],
                 ["ψ_P", [0.0, 1.0]]])

generate_data(X_labels, path; N=N)
# run_PAWN(X_labels, path, run_nr; N=N)