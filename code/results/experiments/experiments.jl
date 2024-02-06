using ArgParse

include("../../model/main.jl")
include("helpers.jl")


global taxtype_dict = Dict(
    :τᴵ => "incometax",
    :τᴷ => "capitaltax",
    :τˢ => "salestax",
    :τᴾ => "profittax",
    :τᴱ => "energytax",
    :τᶜ => "carbontax",
    :prog => "progressivty",
    :green_limit => "green_limit"
)


function getfilepath_tax(
    folderpath::String,
    taxtype::Symbol,
    taxrate::Float64,
    sim_nr::Int64;
    is_agent_data::Bool=false
)
    #folderpath *= "result_data/"
    folderpath = joinpath(pwd(), "OFAT_experiments/result_data/")
    # Convert taxtype Symbol to String name
    folderpath *= taxtype_dict[taxtype]

    if is_agent_data
        folderpath *= "_agentdata"
    end
    #print(folderpath * "_$(taxrate)_$(sim_nr).csv")
    return folderpath * "_$(taxrate)_$(sim_nr).csv"
end

function getfilepath_params(
    folderpath::String,
    param::Symbol,
    param_value::Float64,
    sim_nr::Int64;
    is_agent_data::Bool=false
)
    #folderpath *= "result_data/"
    folderpath = joinpath(pwd(), "OFAT_experiments/result_data/")
    # Convert taxtype Symbol to String name
    folderpath *= param_dict[param]

    if is_agent_data
        folderpath *= "_agentdata"
    end
    #print(folderpath * "_$(taxrate)_$(sim_nr).csv")
    return folderpath * "_$(param_value)_$(sim_nr).csv"
end


# function getfilepath_tax(
#     folderpath::String,
#     taxtype::Symbol,
#     taxrate::Float64,
#     sim_nr::Int64;
#     is_agent_data::Bool=false
# )
#     # Convert taxtype Symbol to String name
#     folderpath = joinpath(folderpath, "result_data")
    
#     folderpath *= taxtype_dict[taxtype]

#     if is_agent_data
#         folderpath *= "_agentdata"
#     end

#     return folderpath * "_$(taxrate)_$(sim_nr).csv"
# end



"""


Runs OFAT experiment on the various tax rates
"""
function OFAT_taxrates(
    folderpath::String,
    taxtypes::Dict;
    modelvars_to_save::Union{Nothing, Vector{Symbol}} = nothing,
    n_per_taxtype::Int64=10,
    n_per_taxrate::Int64=10,
    save_agent_df::Bool=false
)

#if modelvars_to_save is empty list than make it nothing
    if isempty(modelvars_to_save)
        modelvars_to_save = nothing
    end

    for (taxtype, raterange) in taxtypes

        println("type $taxtype, range $raterange")

        for taxrate in Base._linspace(raterange[1], raterange[2], n_per_taxtype)
            
            # Set up array with changed tax rate, to be introduced in period t_warmup
            changed_taxrates = [(taxtype, taxrate)]

            Threads.@threads for sim_nr in 1:n_per_taxrate

                # Set seed according to simulation number
                seed = 1000 + sim_nr

                println("   thread: $(Threads.threadid()), taxrate: $(taxrate), sim nr: $(sim_nr)")

                #Run the model with changed tax rate
                agent_df, model_df = run_simulation(
                    changed_taxrates = changed_taxrates,
                    show_full_output = false,
                    thread_nr = Threads.threadid(),
                    seed = seed,
                    savedata = true
                )
                #replace df with dummy data frames : df = DataFrame(rand(10, 10), :auto)
                # agent_df = DataFrame(rand(10, 10), :auto)
                # model_df = DataFrame(rand(10, 10), :auto)

                # Save agent data

                if save_agent_df
                    outputfilepath_agents = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr, is_agent_data=true)
                    CSV.write(outputfilepath_agents, agent_df)
                end

                # Save model data
                outputfilepath_model = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr)
                if isnothing(modelvars_to_save)
                    # Save full model output
                    CSV.write(outputfilepath_model, model_df)
                else
                    # Save selected model output
                    CSV.write(outputfilepath_model, model_df[!, modelvars_to_save])
                end


                # # Process output data and save to file
                # results = savefulloutput(runoutput, sim_nr; return_as_df=true)
                # outputfilepath = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr)
                # CSV.write(outputfilepath, results)

                # outputfilepath = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr, ishouseholddata=true)
                # CSV.write(outputfilepath, householddata[2])
            end
        end
    end
end


# function getfilepath_gp(
#     folderpath::String,
#     paramtype::Symbol,
#     paramval::Float64,
#     t_introduction::Int64,
#     sim_nr::Int64;
#     isfirmdata::Bool=false,
#     ishouseholddata::Bool=false
# )

#     folderpath *= "$(string(paramtype))_$(paramval)_$(t_introduction)"

#     if isfirmdata
#         folderpath *= "_firmdata"
#     elseif ishouseholddata
#         folderpath *= "_householddata"
#     end

#     return folderpath * "_$(sim_nr).csv"
# end


"""
Runs OFAT experiment on global parameters
"""
function OFAT_globalparam(
    folderpath::String,
    globalparams::Dict;
    modelvars_to_save::Union{Nothing, Vector{Symbol}} = nothing,
    n_per_paramtype::Int64=10,
    n_per_paramval::Int64=10,
    save_agent_df::Bool=false
    )

    #if modelvars_to_save is empty list, then make it nothing
    if isempty(modelvars_to_save)
        modelvars_to_save = nothing
    end

    for (paramtype, paramrange) in globalparams
        println("type $paramtype, range $paramrange")

        for newparamval in Base._linspace(paramrange[1], paramrange[2], n_per_paramtype)
            #newparamval = round(newparamval, digits=2)
            
            # Set up array with changed parameter value, to be introduced in period t_warmup
            changed_params = [(paramtype, newparamval)]

            Threads.@threads for sim_nr in 1:n_per_paramval
                seed = 1000+sim_nr

                println("   thr $(Threads.threadid()), param $(paramtype)=$(newparamval), sim_nr $(sim_nr)")

                # Run the model with changed tax rate
                agent_df, model_df = run_simulation(
                    changed_params=Dict(changed_params),
                    show_full_output=false,
                    thread_nr=Threads.threadid(),
                    seed=seed,
                    savedata=true
                )

                # Save agent data
                if save_agent_df
                    outputfilepath_agents = getfilepath_tax(folderpath,  paramtype, newparamval, sim_nr, is_agent_data=true)
                    CSV.write(outputfilepath_agents, agent_df)
                end

                # Save model data
                outputfilepath_model = getfilepath_tax(folderpath, paramtype, newparamval, sim_nr)
                if isnothing(modelvars_to_save)
                    # Save full model output
                    CSV.write(outputfilepath_model, model_df)
                else
                    # Save selected model output
                    CSV.write(outputfilepath_model, model_df[!, modelvars_to_save])
                end
            end
        end
    end
end


#     Define ranges of tax rates, contains the range and the time step of introduction and duration
    
#     /FIXME Time ranges are not really necessary? 
#     Do we need the definition within the function?
#     globalparams = [
#         # [:p_f, [0.1, 0.3], 420, 12]
#         [:prog, [-1.0, 0.0], 300, Inf]
#     ]

#     for (taxtype, raterange) in taxtypes

#         println("type $taxtype, range $raterange")

#         for taxrate in Base._linspace(raterange[1], raterange[2], n_per_taxtype)
            
#             # Set up array with changed tax rate, to be introduced in period t_warmup
#             changed_taxrates = [(taxtype, taxrate)]

#             Threads.@threads for sim_nr in 1:n_per_taxrate

#                 # Set seed according to simulation number
#                 seed = 1000 + sim_nr

#                 println("   thread: $(Threads.threadid()), taxrate: $(taxrate), sim nr: $(sim_nr)")

#  /FIXME Where do timeranges come from??
#  !!! Was this supposed to do multiple OFATs at once?
#     for (paramtype, paramrange, t_introduction, t_duration) in globalparams

#         println("type $paramtype, range $paramrange")

#         for newparamval in paramrange

#             newparamval = round(newparamval, digits=2)
            
#             # Set up array with changed parameter value, to be introduced in period t_warmup
#             changedparams_ofat = Dict(paramtype => [newparamval, 0., t_introduction, t_duration])

#             Threads.@threads for sim_nr in 1:n_per_paramval

#                 seed = 1000+sim_nr

#                 println("   thr $(Threads.threadid()), param $(paramtype)=$(newparamval), sim_nr $(sim_nr)")

#                 # Run the model with changed tax rate
#                 agent_df, model_df = run_simulation(
#                     changed_params_ofat=changedparams_ofat,
#                     show_full_output=false,
#                     thread_nr=Threads.threadid(),
#                     seed=seed,
#                     savedata=true
#                 )

#                 Save agent data




# function getfilepath_gp_tax(
#     folderpath::String,
#     taxtype::Symbol,
#     taxrate::Float64,
#     paramtype::Symbol,
#     paramval::Float64,
#     t_introduction::Int64,
#     sim_nr::Int64;
#     isfirmdata::Bool=false,
#     ishouseholddata::Bool=false
#     )

#     if taxtype == :τᴵ
#         folderpath *= "incometax"
#     elseif taxtype == :τᴷ
#         folderpath *= "capitaltax"
#     elseif taxtype == :τˢ
#         folderpath *= "salestax"
#     elseif taxtype == :τᴾ
#         folderpath *= "profittax"
#     elseif taxtype == :τᴱ
#         folderpath *= "energytax"
#     else
#         folderpath *= "carbontax"
#     end

#     folderpath *= "_$(taxrate)_$(string(paramtype))_$(paramval)_$(t_introduction)"

#     if isfirmdata
#         folderpath *= "_firmdata"
#     elseif ishouseholddata
#         folderpath *= "_householddata"
#     end

#     return folderpath *= "_$(sim_nr).csv"
# end


# """
# Runs OFAT experiment on global parameters and tax rates simulataneously
# """
# function OFAT_globalparam_taxrates(
#     folderpath::String;
#     n_per_paramtype::Int64=10,
#     n_per_paramval::Int64=10
#     )

#     # Define ranges of tax rates, contains the range and the time step of introduction and duration

#     # taxrates = Dict(
#     #     :τᶜ => (0.3, 0.3)
#     # )

#     globalparams = [
#         [:prog, [-1.0, 0.0], 300, Inf],
#         # [:p_f, [0.3], 420, 12]
#     ]

#     # Define tax rate used for tax experiment
#     taxtype = :τᶜ
#     taxrate = 0.3

#     for (paramtype, paramrange, t_introduction, t_duration) in globalparams

#         for newparamval in paramrange

#             newparamval = round(newparamval, digits=2)
            
#             # Set up array with changed parameter value, to be introduced in period t_warmup
#             changedtaxrates = [(taxtype, taxrate)]
#             changedparams_ofat = Dict(paramtype => [newparamval, 0., t_introduction, t_duration])

#             Threads.@threads for sim_nr in 1:n_per_paramval

#                 seed = 1000+sim_nr

#                 println("   thr $(Threads.threadid()), param $(paramtype)=$(newparamval), sim_nr $(sim_nr)")

#                 # Run the model with changed tax rate
#                 runoutput, firmdata, householddata = run_simulation(
#                     changedparams_ofat=changedparams_ofat,
#                     changedtaxrates=changedtaxrates,
#                     full_output=false,
#                     threadnr=Threads.threadid(),
#                     track_firms_households=true,
#                     seed=seed
#                 )

#                 # Process output data and save to file

#                 results = savefulloutput(runoutput, sim_nr; return_as_df=true)
#                 outputfilepath = getfilepath_gp_tax(folderpath, taxtype, taxrate, paramtype, 
#                                                     newparamval, t_introduction, sim_nr)
#                 CSV.write(outputfilepath, results)

#                 # outputfilepath = getfilepath_gp_tax(folderpath, taxtype, taxrate, paramtype, 
#                 #                                     newparamval, t_introduction, sim_nr, isfirmdata=true)
#                 # CSV.write(outputfilepath, firmdata)

#                 outputfilepath = getfilepath_gp_tax(folderpath, taxtype, taxrate, paramtype, 
#                                                     newparamval, t_introduction, sim_nr, ishouseholddata=true)
#                 CSV.write(outputfilepath, householddata[2])
#             end
#         end
#     end
# end


function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--n_per_taxtype"
            help="number of simulations per tax type"
            arg_type=Int64
            default=20
        "--n_per_taxrate"
            help="number of simulations per tax rate"
            arg_type=Int64
            default=5
        "--outputpath"
            help="path to directory with output files"
            arg_type=String
            default="results/experiments/OFAT_experiments/"
    end

    return parse_args(s)
end


function main()

    # Parse command line arguments
    parsed_args = parse_commandline()
    n_per_type = parsed_args["n_per_taxtype"]
    n_per_val = parsed_args["n_per_taxrate"]
    folderpath = parsed_args["outputpath"]

    # Select tax types and ranges to test
    # taxtype => (τ_start, τ_end, δτ::Int (optional))
    #   τ_start: tax rate from t=0 to t=t_warmup
    #   τ_end: tax rate at t=T
    #   δτ: how many months between increases, default is increase every month (δτ=1).
    # taxrates = Dict(
    #     # :τᴵ => (0.1, 0.6),
    #     # :τᴷ => (0.1, 0.6),
    #     # :τˢ => (0.0, 0.5),
    #     # :τᴾ => (0.1, 0.6),
    #     #:τᴱ => (0.1, 0.8),
    #     :τᶜ => (0.7, 0.71) 
    # )

    paramrange = Dict(
         :prog => (-1.5, 0.5)
    #     :p_f => (0.1, 0.3)
    #     :green_limit => (0.1, 0.11)
    )

    # Select which model variables to save from each simulation
    modelvars_to_save = Symbol[
        # :GDP,               # GDP
        # :U,                 # unemployment rate
        # :em_index,    # carbon emissions index
        # :GINI_W            # Gini coefficient of wealth distribution
        # #add markups
        
    ]

    # Run OFAT experiment for tax rates
    # OFAT_taxrates(
    #     folderpath,
    #     taxrates; 
    #     n_per_taxtype = n_per_type,
    #     n_per_taxrate = n_per_val,
    #     modelvars_to_save = modelvars_to_save
    # )
    OFAT_globalparam(
        folderpath,
        paramrange; 
        n_per_paramtype = n_per_type,
        n_per_paramval = n_per_val,
        modelvars_to_save = modelvars_to_save
    )

    # Run OFAT experiment for global parameters
    # OFAT_globalparam(

    # OFAT_globalparam(folderpath, n_per_paramtype=n_per_type, n_per_paramval=n_per_val)

    # OFAT_globalparam_taxrates(folderpath, n_per_paramtype=n_per_type, n_per_paramval=n_per_val)
end

main()


# #print output of getfolderpath
# folderpath = "jkklkjlkkj"
# taxtype = :τᴱ
# taxrate = 0.1
# sim_nr = 1
# outputfilepath = getfilepath_tax(folderpath, taxtype, taxrate, sim_nr, is_agent_data=true)
# #print(outputfilepath,"///")
# #print(outputfilepath)
# #ommit the last part of the outputfilepath


# #give me the path using pwd() to get current working directory and move into OFAT_experiments folderpath

# #move into OFAT_experiments folder with outputfilepath
# #outputfilepath = joinpath(pwd(), "OFAT_experiments/result_data")



# #add outputfile in csv format to outputfilepath
# #outputfilepath = joinpath(outputfilepath, "outputfile.csv")
# #print current working directory
# print(outputfilepath)
# #create a dummy dataframe
# df = DataFrame(rand(10, 10), :auto)
# #save the dummy dataframe



#CSV.write(outputfilepath, df)
