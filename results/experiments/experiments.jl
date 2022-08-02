using ArgParse

include("../../model/main.jl")
include("experiment_helpers.jl")


function getfilepath(
    folderpath::String,
    taxtype::Symbol,
    taxrate::Float64,
    sim_nr::Int64
    )

    if taxtype == :τᴵ
        folderpath *= "incometax"
    elseif taxtype == :τᴷ
        folderpath *= "capitaltax"
    elseif taxtype == :τˢ
        folderpath *= "salestax"
    elseif taxtype == :τᴾ
        folderpath *= "profittax"
    elseif taxtype == :τᴱ
        folderpath *= "energytax"
    else
        folderpath *= "carbontax"
    end

    return folderpath * "_$(taxrate)_$(sim_nr).csv"
end


"""


Runs OFAT experiment on the various tax rates
"""
function OFAT_taxrates(
    folderpath::String;
    n_per_taxtype::Int64=10,
    n_per_taxrate::Int64=10,
    t_warmup::Int64=100,
    )

    # Define ranges of tax rates
    taxrates = Dict(
        # :τᴵ => (0.0, 0.8),
        # :τᴷ => (0.0, 0.5),
        # :τˢ => (0.0, 0.5),
        # :τᴾ => (0.0, 0.8),
        # :τᴱ => (0.0, 1.0),
        :τᶜ => (0.1, 0.5)
    )

    # lk = Threads.ReentrantLock()

    for (taxtype, raterange) in taxrates

        # firstsave = true

        println("type $taxtype, range $raterange")

        Threads.@threads for taxrate in LinRange(raterange[1], raterange[2], n_per_taxtype)

            taxrate = round(taxrate, digits=2)
            
            # Set up array with changed tax rate, to be introduced in period t_warmup
            changedtaxrates = [(taxtype, taxrate)]

            # results = nothing

            for sim_nr in 1:n_per_taxrate

                seed = 1000+sim_nr

                println("   thr $(Threads.threadid()), tr $(taxrate), sim_nr $(sim_nr)")

                # Run the model with changed tax rate
                runoutput = run_simulation(
                    changedtaxrates=changedtaxrates,
                    full_output=false,
                    threadnr=Threads.threadid(),
                    seed=seed
                )

                # Save results of run
                # if results == nothing
                #     results = savefulloutput(runoutput, sim_nr; return_as_df=true)
                #     results[!, "taxrate"] = [taxrate]
                # else
                #     push!(results, vcat(savefulloutput(runoutput, sim_nr), [taxrate]))
                # end

                # Process output data and save to file
                results = savefulloutput(runoutput, sim_nr; return_as_df=true)
                outputfilepath = getfilepath(folderpath, taxtype, taxrate, sim_nr)
                CSV.write(outputfilepath, results)
            end

            # Save simulation data, wait for other threads already
            # augmenting the file.

            # outputfilepath = getfilepath(folderpath, taxtype, run_nr, seed)
            # CSV.write(outputfilepath, results; append=!firstsave)

            # lock(lk) do
            #     println("$(Threads.threadid()) is saving...")
            #     CSV.write(outputfilepath, results; append=!firstsave)
            #     firstsave = false
            # end
        end
    end
end


function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--n_per_taxtype"
            help="number of simulations per tax type"
            arg_type=Int64
            default=5
        "--n_per_taxrate"
            help="number of simulations per tax rate"
            arg_type=Int64
            default=10
        "--outputpath"
            help="path to directory with output files"
            arg_type=String
            default="results/experiments/OFAT_experiments/"
    end

    return parse_args(s)
end


function main()

    parsed_args = parse_commandline()

    n_per_taxtype = parsed_args["n_per_taxtype"]
    n_per_taxrate = parsed_args["n_per_taxrate"]
    folderpath = parsed_args["outputpath"]

    OFAT_taxrates(folderpath; n_per_taxtype=n_per_taxtype, n_per_taxrate=n_per_taxrate)
end

main()