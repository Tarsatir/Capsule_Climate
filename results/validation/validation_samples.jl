"""
Runs MC simulation to acquire data for validation excercises.
"""

using PyCall
using ArgParse

include("../../model/main.jl")
include("../../results/experiments/experiment_helpers.jl")


function runMCreplication(
    repl_i::Int64;
    t_warmup::Int64=300
    )

    println("replication number $repl_i")

    # Run simulation with default (calibrated) parameters
    runoutput = run_simulation(
        full_output=false
    )

    # Save output data
    df = DataFrame(
        :GDP => runoutput.GDP[t_warmup:end],
        :GDP_growth => runoutput.GDP_growth[t_warmup:end],
        :C => runoutput.C[t_warmup:end],
        :I => runoutput.I[t_warmup:end],
        :wages => runoutput.wages[t_warmup:end],
        :prices => runoutput.prices[t_warmup:end],
        :markups => runoutput.markups[t_warmup:end],
        :TotDebt => runoutput.TotDebt[t_warmup:end],
        :EnDem => runoutput.EnDem[t_warmup:end],
        :U => runoutput.U[t_warmup:end],
        :LIS => runoutput.LIS[t_warmup:end],
        :Em => runoutput.emissions_total[t_warmup:end],
        :EmIndex => runoutput.emissions_index[t_warmup:end],
        :RD => runoutput.RD[t_warmup:end],
        :inventories => runoutput.inventories[t_warmup:end]
    )

    filepath = "results/validation/validation_samples/valoutput_$repl_i.csv"

    CSV.write(filepath, df)
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
            default=10
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

    for repl_i in start:finish
        runMCreplication(repl_i)
    end
end

main()