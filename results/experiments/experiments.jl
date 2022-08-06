using ArgParse

include("../../model/main.jl")
include("experiment_helpers.jl")


function getfilepath_tax(
    folderpath::String,
    taxtype::Symbol,
    taxrate::Float64,
    sim_nr::Int64;
    isfirmdata::Bool=false,
    ishouseholddata::Bool=false
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

    if isfirmdata
        folderpath *= "_firmdata"
    elseif ishouseholddata
        folderpath *= "_householddata"
    end

    return folderpath * "_$(taxrate)_$(sim_nr).csv"
end


"""


Runs OFAT experiment on the various tax rates
"""
function OFAT_taxrates(
    folderpath::String;
    n_per_taxtype::Int64=10,
    n_per_taxrate::Int64=10
    )

    # Define ranges of tax rates
    taxrates = Dict(
        # :τᴵ => (0.0, 0.8),
        # :τᴷ => (0.0, 0.5),
        # :τˢ => (0.0, 0.5),
        # :τᴾ => (0.0, 0.8),
        :τᴱ => (0.1, 0.7),
        # :τᶜ => (0.1, 0.5)
    )

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
                runoutput, firmdata, householddata = run_simulation(
                    changedtaxrates=changedtaxrates,
                    full_output=false,
                    threadnr=Threads.threadid(),
                    track_firms_households=true,
                    seed=seed
                )

                # Process output data and save to file
                results = savefulloutput(runoutput, sim_nr; return_as_df=true)
                outputfilepath = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr)
                CSV.write(outputfilepath, results)

                outputfilepath = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr, isfirmdata=true)
                CSV.write(outputfilepath, firmdata)

                outputfilepath = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr, ishouseholddata=true)
                CSV.write(outputfilepath, householddata[2])
            end
        end
    end
end


function getfilepath_gp(
    folderpath::String,
    paramtype::Symbol,
    paramval::Float64,
    t_introduction::Int64,
    sim_nr::Int64;
    isfirmdata::Bool=false,
    ishouseholddata::Bool=false
    )

    folderpath *= "$(string(paramtype))_$(paramval)_$(t_introduction)"

    if isfirmdata
        folderpath *= "_firmdata"
    elseif ishouseholddata
        folderpath *= "_householddata"
    end

    return folderpath * "_$(sim_nr).csv"
end


"""
Runs OFAT experiment on global parameters
"""
function OFAT_globalparam(
    folderpath::String;
    n_per_paramtype::Int64=10,
    n_per_paramval::Int64=10
    )

    # Define ranges of tax rates, contains the range and the time step of introduction and duration
    globalparams = [
        [:p_f, [0.1, 0.3], 420, 12]
    ]

    for (paramtype, paramrange, t_introduction, t_duration) in globalparams

        println("type $paramtype, range $paramrange")

        Threads.@threads for newparamval in repeat(paramrange, n_per_paramtype)

            newparamval = round(newparamval, digits=2)
            
            # Set up array with changed parameter value, to be introduced in period t_warmup
            changedparams_ofat = Dict(paramtype => [newparamval, 0., t_introduction, t_duration])

            for sim_nr in 1:n_per_paramval

                seed = 1000+sim_nr

                println("   thr $(Threads.threadid()), param $(paramtype)=$(newparamval), sim_nr $(sim_nr)")

                # Run the model with changed tax rate
                runoutput, firmdata, householddata = run_simulation(
                    changedparams_ofat=changedparams_ofat,
                    full_output=false,
                    threadnr=Threads.threadid(),
                    track_firms_households=true,
                    seed=seed
                )

                # Process output data and save to file
                results = savefulloutput(runoutput, sim_nr; return_as_df=true)
                outputfilepath = getfilepath_gp(folderpath, paramtype, newparamval, t_introduction, sim_nr)
                CSV.write(outputfilepath, results)

                outputfilepath = getfilepath_gp(folderpath, paramtype, newparamval, t_introduction, sim_nr, isfirmdata=true)
                CSV.write(outputfilepath, firmdata)

                outputfilepath = getfilepath_gp(folderpath, paramtype, newparamval, t_introduction, sim_nr, ishouseholddata=true)
                CSV.write(outputfilepath, householddata[2])
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
            default=5
        "--n_per_taxrate"
            help="number of simulations per tax rate"
            arg_type=Int64
            default=30
        "--outputpath"
            help="path to directory with output files"
            arg_type=String
            default="results/experiments/OFAT_experiments/"
    end

    return parse_args(s)
end


function main()

    parsed_args = parse_commandline()

    n_per_type = parsed_args["n_per_taxtype"]
    n_per_val = parsed_args["n_per_taxrate"]
    folderpath = parsed_args["outputpath"]

    OFAT_taxrates(folderpath; n_per_taxtype=n_per_type, n_per_taxrate=n_per_val)

    # OFAT_globalparam(folderpath, n_per_paramtype=n_per_type, n_per_paramval=n_per_val)
end

main()