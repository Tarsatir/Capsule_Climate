using PyCall

include("model/main.jl")

seed = 2345

# Ensure the output directory exists
output_dir = "./within_time_ofat"
isdir(output_dir) || mkdir(output_dir)

# Run model n times with different seeds and save data in "within_time_ofat" folder
for i in 1:9
    local_seed = (i == 1) ? seed : seed + i
    
    run_simulation(
        T = 660;
        savedata = true,
        show_full_output = true,
        showprogress = true,
        seed = local_seed,
        changed_taxrates = [(:τᶜ, 0.4)]
    )
    
    name = "$(local_seed)_model.csv"
    
    cp(name, joinpath(output_dir, name), force=true)
    
    # Delete the original file
    rm(name)
end





nothing
