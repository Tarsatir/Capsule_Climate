using PyCall
using DataFrames
using LinearAlgebra

include("../../model/main.jl")
include("../../results/experiments/experiment_helpers.jl")

spopt = pyimport("scipy.optimize")


function run_model(
    xcurr::Vector{Float64}, 
    labels::Vector{String},
    truemoments::DataFrame,
    W::Matrix{Float64},
    rescale::Vector{Float64}
    )

    changedparams = Dict(labels .=> xcurr)
    runoutput = run_simulation(
        changed_params=changedparams,
        full_output=true
    )

    truemoments_values = collect(values(truemoments[1, :]))
    # println(truemoments_values)

    sim_moments = convertrunoutput(runoutput, 0; return_as_df=true)
    simmoments_values = collect(values(sim_moments[1, propertynames(truemoments)]))

    # println(simmoments_values)
    
    crit = compute_criterion(
                simmoments_values,
                truemoments_values,
                W,
                rescale
            )

    return crit
end

function compute_criterion(
    simmoments::Vector{Float64},
    truemoments::Vector{Float64},
    W::Matrix{Float64},
    rescale::Vector{Float64}
    )

    # println(simmoments)
    # println(truemoments)

    simmoments .*= rescale
    truemoments .*= rescale

    error = simmoments .- truemoments
    return transpose(error) * W * error
end


function min_MSE(xlabels, xbounds, truemoments, rescale)

    options = Dict("maxfev"=> 100, "disp"=> true)

    W = Matrix{Float64}(I, length(rescale), length(rescale))

    optimimal_θ = spopt.minimize(
        run_model,
        x0 = [0.8, 0., 0.8, 0.7, 0.005, 0.2],
        args = Tuple([xlabels, truemoments, W, rescale]),
        method = "Nelder-Mead",
        bounds = xbounds,
        tol = 1e-8,
        options = options
    )

    df = DataFrame(
        :errors => optimimal_θ["final_simplex"][end]
    )

    println(optimimal_θ)
    println(optimimal_θ["x"])

    CSV.write("parameters/calibration/opt_errors.csv", df)
end


function run_MSM_calibration()

    xlabels = ["α_cp", "prog", "ω", "λ", "κ_upper", "p_f"]
    xbounds = [[0.6, 1.0],
               [-1.0, 1.0],
               [0.0, 1.0],
               [0.0, 1.0],
               [0.0, 0.05],
               [0.0, 1.0]]

    truemoments = DataFrame(CSV.File("parameters/calibration/true_moments.csv"))

    rescale = Float64[100., 10., 0., 1., 100., 0., 0., 0., 10., 1., 1., 1.]

    min_MSE(xlabels, xbounds, truemoments, rescale)
end

run_MSM_calibration()