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
function generate_data()

    # Set up which parameter to apply GSA to, and the range of the parameter
    X_labels = Dict([
        ["α_cp", [0.6, 1.0]],
        ["μ1", [0.0, 0.5]]
        ])

    # Number of uncertain parameters
    M = length(X_labels)

    # println(X_labels)

    # Define parameter distributions
    samp_strat = "lhs"
    N = 100

    # Call SAMP function to get the input parameters for running the model
    py"""

    def call_AAT_sampling(samp_strat, M, X_labels, N):
        
        import scipy.stats as stats
        import numpy as np
        from SAFEpython.sampling import AAT_sampling
        
        # Define distribution of parameters
        distr_fun = [stats.uniform] * M

        distr_par = [np.nan] * M
        for i,key in enumerate(X_labels):
            distr_par[i] = [X_labels[key][0], X_labels[key][1] - X_labels[key][0]]

        X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
        
        return X

    """

    X = py"call_AAT_sampling"(samp_strat, M, X_labels, N)
    # display(X)
    
    # Run simulations for the given parameter values
    Y = zeros(N)

    for i in 1:N

        println("Simulation nr: $i")

        # Retrieve variables, repackage to pass to simulation
        Xi = X[i,:]
        changed_params = Dict(k=>Xi[i] for (i,k) in enumerate(keys(X_labels)))

        # println(changed_params)

        # Run simulation
        Yi = run_simulation(
                changed_params=changed_params,
                full_output=false
            )

        # Store output data
        Y[i] = Yi
    end

    # Set up dataframe containing results
    df = DataFrame(Dict(x=>X[:,i] for (i,x) in enumerate(keys(X_labels))))
    df[!, "GDP"] = Y

    # Write results to csv
    CSV.write("parameters/sensitivity/sensitivity_runs/sensitivity_run_1.csv", df)

end


"""
Calls SAFE toolbox in Python script.
"""
function run_PAWN()

    # Include Python file containing GSA functions
    @pyinclude("parameters/sensitivity/run_GSA.py")

    # Call PAWN function
    py"run_PAWN"()

end

generate_data()
# run_PAWN()