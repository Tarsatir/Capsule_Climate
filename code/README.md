### README
This is the collection of code needed to replicate all results in the paper "Carbon Pricing Drives Critical Transition to Green Growth" by Isaak Mengesha Students under his supervision at the University of Amsterdam. The code is written in `Julia` and `Python`. The `Julia` code is used to run the model, and the `Python` code is used to produce the plots.

### Model description
The function 'run_simulation()' called in "run_model.jl" allows for the specification of certain number of experiments. If you want to merely toy with a single model run, you can do so by adjusting the type of parameters passed to the function, and when and how they are changed. The function 'run_simulation()' is specified in main.jl. For certain fixed parameter changs make adjustments to the global_parameter.jl (note that not all relevant parameters can be found there, as some are agent specific). 

The sensitivity analysis for the model was done with the 'run_sensitivity_analysis()' function. This function is specified in sensitivity_analysis.jl. For a fixed number of parameters and corresponding ranges the PAWN algorithm is implemented. 

### Running the model

Run the model by starting a `Julia` session in the main folder, and typing the following command:

```julia
>julia include("run_model.jl")
```
This command will run the model for the default settings, and produce the output plots.
In order to replicate the main finding of the paper (the critical transition), you will need to run the model for a large number of times for varying 'green_limits' using the 'changed_params' option in the 'run_simulation()' function. 

### Plots and Results
For some tentative plots execute 'plot_macro_vars.py' these will not replicate the plots in the paper exactly, as the data generated for the paper is not included in this repository.