using ArgParse

include("../../model/main.jl")
include("experiment_helpers.jl")


function getfilepath(
    folderpath::String,
    taxtype::Symbol
    )

    if taxtype == :τᴵ
        return folderpath * "incometax.csv"
    elseif taxtype == :τᴷ
        return folderpath * "capitaltax.csv"
    elseif taxtype == :τˢ
        return folderpath * "salestax.csv"
    elseif taxtype == :τᴾ
        return folderpath * "profittax.csv"
    elseif taxtype == :τᴱ
        return folderpath * "energytax.csv"
    else
        return folderpath * "carbontax.csv"
    end
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
        :τᶜ => (0.0, 0.2)
    )

    lk = Threads.ReentrantLock()

    for (i, (taxtype, raterange)) in enumerate(taxrates)

        firstsave = true
        outputfilepath = getfilepath(folderpath, taxtype)

        println("type $taxtype, range $raterange")

        Threads.@threads for taxrate in LinRange(raterange[1], raterange[2], n_per_taxtype)
            changedtaxrates = [(taxtype, taxrate)]

            results = nothing

            for i in 1:n_per_taxrate

                println("   thr $(Threads.threadid()), tr $(taxrate), i $(i)")

                # Run the model with changed tax rate
                runoutput = run_simulation(
                    changedtaxrates=changedtaxrates,
                    full_output=false,
                    threadnr=Threads.threadid()
                )

                # Save results of run
                if results == nothing
                    results = convertrunoutput(runoutput; return_as_df=true)
                    results[!, "taxrate"] = [taxrate]
                else
                    push!(results, vcat(convertrunoutput(runoutput), [taxrate]))
                end
            end

            # Save simulation data, wait for other threads already
            # augmenting the file.
            lock(lk) do
                println("$(Threads.threadid()) is saving...")
                CSV.write(outputfilepath, results; append=!firstsave)
                firstsave = false
            end
        end
    end
end


function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin
        
        "--n_per_taxtype"
            help="number of simulations per tax type"
            arg_type=Int64
            default=10
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