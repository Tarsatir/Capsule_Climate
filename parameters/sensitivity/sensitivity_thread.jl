include("../../model/main.jl")

path(run_nr, pathroot, thread_nr) = "$(pathroot)sensitivity_run_$(run_nr)_thr_$(thread_nr).csv"

"""
    function runsensitivitythread(;n_per_epoch::Int = 10)

Runs the sensitivity thread based on the number of thread number, given 
    by comand line argument. Code changed to be run in parallel on the HPC.
"""
function runsensitivitythread(;
    n_per_epoch::Int = 3,
    run_nr::Int = 5,
    M::Int = 9
)

    # Define path of used data file
    pathroot = ARGS[1]
    threadnr = parse(Int64, ARGS[2])
    threadpath = path(run_nr, pathroot, threadnr)

    df = DataFrame(CSV.File(threadpath))
    n_per_thread = nrow(df)

    # Determine total amount of epochs
    n_epochs = ceil(Int64, n_per_thread / n_per_epoch)

    # Allow n epochs to be changed for testing purposes
    if length(ARGS) > 2
        n_epochs = min(parse(Int64, ARGS[3]), n_epochs)
    end

    total_completed_runs = 0
    
    for _ in 1:n_epochs
        df = DataFrame(CSV.File(threadpath))

        firstrow = total_completed_runs + 1
        lastrow = firstrow + n_per_epoch - 1

        labels = names(df)[1:M]

        for row in firstrow:lastrow

            # Generate dictionary with changed parameters
            changedparams = Dict()
            for label in labels
                changedparams[label] = df[row, label]
            end

            # Run the model with changed parameters
            GDP_g, GINI_I, GINI_W, U, FGT = run_simulation(
                changed_params=changedparams,
                full_output=false;
                threadnr=threadnr
            )

            # Write GDP data to dataframe
            df[row, "GDP_1st"] = mean(GDP_g)
            df[row, "GDP_2nd"] = var(GDP_g)
            df[row, "GDP_3rd"] = skewness(GDP_g)
            df[row, "GDP_4th"] = kurtosis(GDP_g)

            # Write Gini data to dataframe
            df[row, "Gini_I_1st"] = mean(GINI_I)
            df[row, "Gini_I_2nd"] = var(GINI_I)

            df[row, "Gini_W_1st"] = mean(GINI_W)
            df[row, "Gini_W_2nd"] = var(GINI_W)

            # Write unemployment data to dataframe
            df[row, "U_1st"] = mean(U)
            df[row, "U_2nd"] = var(U)

            # Write poverty data to dataframe
            df[row, "FGT_1st"] = mean(FGT)
            df[row, "FGT_2nd"] = var(FGT)

            # Update total completed runs
            total_completed_runs += 1
            if total_completed_runs == nrow(df)
                break
            end

        end

        # At end of epoch, write (intermediate) results to dataframe
        CSV.write(threadpath, df)

        if total_completed_runs == nrow(df)
            break
        end
    end

end

runsensitivitythread()