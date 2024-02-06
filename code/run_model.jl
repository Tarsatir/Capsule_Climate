using PyCall

include("model/main.jl")

seed = 1234

run_simulation(
    T = 660;
    savedata = true,
    show_full_output = true,
    showprogress = true,
    seed = 1234,
    save_firmdata = true
)

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
#     changed_params=Dict([(:p_f, 0.2)]),
# )

nothing

# @pyinclude("plotting/plot_macro_vars.py")