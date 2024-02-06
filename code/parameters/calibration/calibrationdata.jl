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
        return folderpath * "cal_output_run$(run_nr)_thread$(thread_nr).csv"
    end
end


function gen_calibration_data(
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
                res = savefulloutput(runoutput, sim_nr; return_as_df=true)
            else
                push!(res, (savefulloutput(runoutput, sim_nr)))
            end

            # If end of epoch is reached, write results to output csv
            if nrow(res) == n_per_epoch || i == n_per_thread
                # println("tread $(thread_nr) is saving...")
                CSV.write(outputfilepath, res; append=i≠n_per_epoch)
                res = nothing
            end

            if i == n_per_thread
                break
            end

        end
    end
end


function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin
        
        "--n_sims", "-N"
            help="number of simulations"
            arg_type=Int64
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
            default="parameters/calibration/calibrationdata/"
    end

    return parse_args(s)
end


function main()

    run_nr = 9

    X_labels = Dict([
        ["α_cp", [0.4, 1.0]],
        ["prog", [-1.0, 1.0]],
        # ["μ1", [0.0, 0.5]],
        ["ω", [0.0, 1.0]],
        ["ϵ", [0.0, 0.1]],
        ["κ_upper", [0.0, 0.05]],
        # ["ψ_E", [0.0, 0.1]],
        ["p_f", [0.0, 1.0]]]
   )

    parsed_args = parse_commandline()

    # Unpack parsed arguments
    n_sims = parsed_args["n_sims"]
    n_per_epoch = parsed_args["n_per_epoch"]
    inputpath = parsed_args["inputpath"]
    outputpath = parsed_args["outputpath"]

    # Specify number of threads and n per thread
    n_threads = Threads.nthreads()
    n_per_thread = ceil(Int64, n_sims / n_threads)

    gen_calibration_data(
        X_labels,
        n_threads,
        n_per_epoch,
        n_per_thread,
        inputpath,
        outputpath,
        run_nr
    )
end

main()