### README
This is the collection of code needed to replicate all results in the paper "Carbon Pricing Drives Critical Transition to Green Growth" by Isaak Mengesha et al. The modeling framework and much of the Julia codebase relies on work by J. Akkerman, developed as part of their Master’s thesis at the University of Amsterdam’s Computational Science Lab (Aug. 2022). Supervised by Isaak Mengesha MSc, and assessed by Dr. ir. Florian Wagener and Dr. Debraj Roy, this thesis provides the methodological cornerstone for the simulations in this repository.

The codebase is written in Julia and Python. The Julia code is used to run the model, and the Python code is used to produce the plots.

### Model description
The function 'run_simulation()' called in "run_model.jl" allows for the specification of certain number of experiments. If you want to merely toy with a single model run, you can do so by adjusting the type of parameters passed to the function, and when and how they are changed. The function 'run_simulation()' is specified in main.jl. For certain fixed parameter changs make adjustments to the global_parameter.jl (note that not all relevant parameters can be found there, as some are agent specific). 

The sensitivity analysis for the model was done with the 'run_sensitivity_analysis()' function. This function is specified in sensitivity_analysis.jl. For a fixed number of parameters and corresponding ranges the PAWN algorithm is implemented. 

### Reproduce results

Check if you have the necessary requirements 'requirements.txt' by running the following command: 'pip-check-reqs'.
The python script 'run_pipeline' executes all relevant julia and python functions to run the experiments and produce the plots in the paper. 
Should you want to run your own experiments navigate to the 'code' folder and make adjustments to the 'run_model.jl' file.

