using CSV
using DataFrames

path = "sensitivity_run_6_thr_1.csv"
output_path = "output_thread1.csv"
# df_input = DataFrame(CSV.File(path))

function writeread(
    n_per_epoch = 50
    )

    res = nothing

    for (i,row) in enumerate(CSV.Rows(path))

        # println(i)

        if (i - 1) % n_per_epoch == 0
            # println("yeet")
            res = DataFrame(
                :μ1 => row.μ1,
                :ϵ => row.ϵ
            )
        else
            push!(res, [row.μ1, row.ϵ])
        end

        if nrow(res) == n_per_epoch 
            CSV.write(output_path, res; append=i≠1)
        end
    end
end

@time writeread()