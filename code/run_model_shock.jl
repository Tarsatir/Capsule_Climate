using PyCall
using FileIO

function move_and_rename_file(source_folder::String, old_filename::String, new_filename::String)
    source_path = joinpath(source_folder, old_filename)
    destination_path = joinpath(source_folder, new_filename)
    
    try
        mv(source_path, destination_path, force=true)
        println("File moved and renamed successfully.")
    catch err
        println("Error: $err")
    end
end

function move_and_rename(old_path::String, new_path::String)
    # Move and rename the folder
    mv(old_path, new_path, force=true)

    # Create a new folder with the old name
    mkdir(old_path)
end








include("model/main.jl")

# run_simulation(
#     T = 660;
#     savedata = true,
#     show_full_output = true,
#     showprogress = true,
#     seed = 1234
# )

# run_simulation(
#     T = 660;
#     savedata = true,
#     show_full_output = true,
#     showprogress = true,
#     seed = seed,
#     changed_taxrates = [(:τᶜ, 0.8)]
# )

# run_simulation(
#     T = 660;
#     savedata = true,
#     show_full_output = true,
#     showprogress = true,
#     seed = seed,
#     changed_taxrates = [(:τᶜ, 0.8, 0.85)],
#     save_firmdata = true
# )
seed = 1234
run_simulation(
    T = 660;
    savedata = true,
    show_full_output = true,
    showprogress = true,
    seed = seed,
    changed_params=Dict([(:p_f, 0.1)]),
)

# Example usage:
source_folder = joinpath(pwd(), "data")
old_filename = "1234_model.csv"
new_filename = "priceshocks/down_priceshock.csv"

# Check if file exists
if isfile(joinpath(source_folder, old_filename))
    println("File exists.")
end

move_and_rename_file(source_folder, old_filename, new_filename)


# Example usage
old_folder = "household_data.csv"
new_folder = "priceshocks/household_data_down.csv"

move_and_rename_file(source_folder, old_folder, new_folder)


run_simulation(
    T = 660;
    savedata = true,
    show_full_output = true,
    showprogress = true,
    seed = seed,
    changed_params=Dict([(:p_f, 0.3)]),
)
source_folder = joinpath(pwd(), "data")
old_filename = "1234_model.csv"
new_filename = "priceshocks/up_priceshock.csv"

move_and_rename_file(source_folder, old_filename, new_filename)

# Example usage:
old_folder = "household_data.csv"
new_folder = "priceshocks/household_data_up.csv"

move_and_rename_file(source_folder, old_folder, new_folder)


run_simulation(
    T = 660;
    savedata = true,
    show_full_output = true,
    showprogress = true,
    seed = seed,
    changed_params=Dict([(:p_f, 0.2)]),
)
source_folder = joinpath(pwd(), "data")
old_filename = "1234_model.csv"
new_filename = "priceshocks/no_priceshock.csv"

move_and_rename_file(source_folder, old_filename, new_filename)

# Example usage:
old_folder = "household_data.csv"
new_folder = "priceshocks/household_data.csv"

move_and_rename_file(source_folder, old_folder, new_folder)

nothing

# @pyinclude("plotting/plot_macro_vars.py")