## Sensitivity analysis

```
julia parameters/sensitivity/sensitivity.jl --n_sims --n_per_epoch --path
```

```
julia parameters/sensitivity/sensitivity.jl --help
```

Alternatively, it may be faster to start a Julia REPL and run the code from here. This can be done in the following way, first type `julia` in the terminal and then once the REPL has started, run the following command:

```julia
julia> include("parameters/sensitivity/sensitivity.jl --n_sims --n_per_epoch --path")
```