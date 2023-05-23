using PyCall

include("model/main.jl")

seed = 1234

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
#     changed_taxrates = [(:τᶜ, 0.5)]
# )

run_simulation(
    T = 1600;
    savedata = true,
    show_full_output = true,
    showprogress = true,
    seed = seed,
    changed_taxrates = [(:τᶜ, 0.2, 0.8)]
)


nothing

# @pyinclude("plotting/plot_macro_vars.py")