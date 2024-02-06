"""
Runs MC simulation to acquire data for validation excercises.
"""

using PyCall
using ArgParse

include("../../model/main.jl")
include("../../results/experiments/helpers.jl")


function runMCreplication(
    repl_i::Int64
    # t_warmup::Int64=300
    )

    println("replication number $repl_i")

    seed = repl_i + 1000

    # Run simulation with default (calibrated) parameters
    runoutput, firmdata, householddata = run_simulation(
        full_output=false,
        seed=seed,
        track_firms_households=true
    )

    df = savefulloutput(runoutput, repl_i; return_as_df=true)

    filepath_runoutput = "results/validation/validation_samples/valoutput_$(repl_i).csv"
    CSV.write(filepath_runoutput, df)

    filepath_firmdata = "results/validation/validation_samples/valfirmdata_$(repl_i).csv"
    CSV.write(filepath_firmdata, firmdata)

    filepath_householddata = "results/validation/validation_samples/householddata_$(repl_i).csv"
    CSV.write(filepath_householddata, householddata[2])
end


function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--proc_i"
            help="index of process"
            arg_type=Int64
            default=1
        "--n_proc"
            help="number of processes"
            arg_type=Int64
            default=1
        "--n_repl"
            help="number of replications"
            arg_type=Int64
            default=100
    end

    return parse_args(s)
end


function main()

    parsed_args = parse_commandline()

    proc_i = parsed_args["proc_i"]
    n_proc = parsed_args["n_proc"]
    n_repl = parsed_args["n_repl"]

    n_per_proc = ceil(Int64, n_repl / n_proc)
    start = (proc_i - 1) * n_per_proc + 1
    finish = (proc_i) * n_per_proc

    Threads.@threads for repl_i in start:finish
        runMCreplication(repl_i)
    end
end

main()