using PyCall
using ArgParse
using DataFrames
using CSV
using Random
using Tables

include("../../model/main.jl")

pysm = pyimport("statsmodels.api")


function run_calibration(
    ytrue::Matrix{Float64},
    outputpath::String
    )

    labels = ["κ_upper", "ω", "ϵ", "α_cp", "p_f", "prog"]
    xbounds = [[0.001, 0.01], [0.0, 1.0], [0.05, 0.5], [0.65, 0.95], [0.1, 1.0], [-1.0, 1.0]]
    
    # Define initial x and weights used for transformations
    xcurr = Float64[0.005227455, 0.429497195, 0.085293161, 0.808659944, 0.340361401, 0.116886424]
    # xcurr = Float64[0.005, 0.5, 0.2, 0.85, 0.3, 0.0]
    xweights = Float64[0.001, 0.1, 0.001, 0.1, 0.01, 0.01]

    # Get y and likelihood from initial x
    y = run_model(xcurr, labels)
    lhcurr = compute_likelihood(y, ytrue)

    filepath = outputpath * "bayes_sampled_params_thread.csv"
    savedata(0, lhcurr, xcurr, labels, filepath; append=false)

    print(lhcurr)

    for i in 1:100

        # Draw new set of parameter values x
        xnew = gen_x(xcurr, xweights, xbounds)
        ynew = run_model(xnew, labels)
        lhnew = compute_likelihood(ynew, ytrue)

        # Compute acceptance chance
        α = min(1, lhnew / lhcurr)

        println(xnew)
        println(lhnew)
        println(α)

        if Random.rand() < α
            # Save parameter set
            xcurr .= xnew
            lhcurr = lhnew
            savedata(i, lhcurr, xcurr, labels, filepath)
            # res = vcat([i, lhcurr], xcurr)
            # CSV.write(filepath, Tables.table(res))
        end
    end
end

function savedata(
    i::Int64,
    lhcurr::Float64,
    xcurr::Vector{Float64},
    labels::Vector{String},
    filepath::String;
    append::Bool=true
    )

    # Initialize dataframe that holds parameter values, write to CSV
    df = DataFrame(vcat(["sim_nr", "lh"], labels) .=> vcat([i, lhcurr], xcurr))
    CSV.write(filepath, df, delim=';', append=append)
end


function compute_likelihood(
    y::Matrix{Float64},
    ytrue::Matrix{Float64}
    )

    f = pysm.nonparametric.KDEMultivariate(
        data=y,
        var_type="ccc",
        bw="normal_reference"
    )
    likelihood = prod(f.pdf(ytrue))
    return likelihood
end


function gen_x(xcurr, xweights, xbounds)

    # Draw which parameter will be changed
    choice = zeros(length(xcurr))
    choice[rand(1:length(xcurr))] = 1

    # Draw stochastic shocks
    ε = Random.randn(length(xcurr)) .* xweights .* choice
    xnew = xcurr .+ ε

    # Make sure parameter values stay within their bounds
    while !inbounds(xnew, xbounds)
        ε = Random.randn(length(xcurr)) .* xweights .* choice
        xnew = xcurr .+ ε
    end
    return xnew
end


function inbounds(
    xnew, 
    xbounds
    )

    for (x, bound) in zip(xnew, xbounds)
        if x < bound[1] || x > bound[2]
            return false
        end
    end
    return true
end


function run_model(
    xcurr::Vector{Float64}, 
    labels::Vector{String}
    )

    changedparams = Dict(labels .=> xcurr)

    runoutput = run_simulation(
        changed_params=changedparams,
        full_output=false;
        threadnr=Threads.threadid()
    )

    return hcat(
            runoutput.GDP_growth[100:end], 
            100 .* runoutput.U[100:end], 
            runoutput.emissions_index[100:end]
           )
end

function main()

    inputpath = "../sensitivity/sensitivity_runs/input_data/"
    outputpath = "calibrationdata/"
    n_threads = 16
    run_nr = 9

    ytrue = Matrix(DataFrame(CSV.File(outputpath * "truedata.csv")))

    run_calibration(ytrue, outputpath)

end

main()