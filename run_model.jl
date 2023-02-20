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
    T = 660;
    savedata = true,
    show_full_output = true,
    showprogress = true,
    seed = seed,
    changed_taxrates = [(:τᶜ, 0.1, 0.5)]
)

nothing

# @pyinclude("plotting/plot_macro_vars.py")