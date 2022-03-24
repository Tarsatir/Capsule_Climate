mutable struct GlobalParam
    # Determine technical innovation process
    ν :: Float64                    # R&D inv propensity
    ξ :: Float64                    # R&D allocation to IN
    ζ :: Float64                    # firm search capabilities
    α1 :: Float64                   # 1st beta dist param for IN
    β1 :: Float64                   # 2nd beta dist param for IN
    κ_lower :: Float64              # 1st beta dist support
    κ_upper :: Float64              # 2nd beta dist support


    γ :: Float64                    # new custommer sample parameter
    μ1 :: Float64                   # kp markup rule
    r :: Float64                    # interest rate
    ι :: Float64                    # desired inventories
    b :: Int                        # payback period
    η :: Int                        # physical scrapping age

    # υ :: Float64                    # markup coefficient
    # ω :: Float64                    # competetiveness weights
    # χ :: Float64                    # replicator dynamics coeff

    Λ :: Float64                    # max debt/sales ratio, regular debt
    # Λ_max :: Float64                # max debt/sales ratio, including credit

    # Determine entrant composition
    φ1 :: Float64                   # 1st Uniform dist support, cp entrant cap
    φ2 :: Float64                   # 2nd Uniform dist support, cp entrant cap
    φ3 :: Float64                   # 1st Uniform dist support, cp entrant liq
    φ4 :: Float64                   # 2nd Uniform dist support, cp entrant liq
    α2 :: Float64                   # 1st beta dist param for kp entrant
    β2 :: Float64                   # 2nd beta dist param for kp entrant

    cu :: Float64                   # capacity utilization for cp
    ϵ :: Float64                    # minimum desired wage increase rate
    max_g_wᴼ :: Float64             # max growth rate of offered wages
    Kg_max :: Float64               # maximum capital growth rate

    # Determine expectation updating cp
    ωD :: Float64                   # memory parameter cp demand estimation
    ωQ :: Float64                   # memory parameter cp quantity estimation
    ωL :: Float64                   # memory parameter cp labor supply estimation
    ωW :: Float64                   # memory parameter cp offered wage estimation

    # Determine household consumption
    α_cp :: Float64                 # parameter controlling MPC of consumers
    c_L_max :: Float64              # maximum share consumed on luxury goods
    a_σ :: Float64                  # 1st parameter governing logistic function
    b_σ :: Float64                  # 2nd parameter governing logistic function

    # Determine household switching
    ψ_E :: Float64                  # chance of employed worker looking for a better paying job
    ψ_Q :: Float64                  # chance of household switching away from cp when demand constrained
    ψ_P :: Float64                  # chance of household switching to cp with better price
    # Nᵈ_share::Float64
end


function initialize_global_params(
        changed_params
    )

    global_param = GlobalParam(
        0.04,                       # ν: R&D inv propensity
        0.5,                        # ξ: R&D allocation to IN
        0.3,                        # ζ: firm search capabilities
        3.0,                        # α1: 1st beta dist param for IN
        3.0,                        # β1: 2nd beta dist param for IN
        -0.03,                      # κ_lower: 1st beta dist support
        0.03,                       # κ_upper: 2nd beta dist support

        0.5,                        # γ: new custommer sample parameter
        0.15,                       # μ1: kp markup rule
        0.0,                        # r: (annual) interest rate
        0.1,                        # ι: desired inventories
        3,                          # b: payback period
        20,                         # η: physical scrapping age

        # 0.05,                       # υ: markup coefficient
        # 1.0,                        # ω: competetiveness weights
        # 1.0,                        # χ: replicator dynamics coeff

        2.0,                        # Λ: max debt/sales ratio regular debt

        # 2.5,                        # Λ_max: max debt/sales ratio including credit

        0.1,                        # φ1: 1st Uniform dist support, cp entrant cap
        0.9,                        # φ2: 2nd Uniform dist support, cp entrant cap
        0.1,                        # φ3: 1st Uniform dist support, cp entrant liq
        0.9,                        # φ4: 2nd Uniform dist support, cp entrant liq
        2.0,                        # α2: 1st beta dist param for kp entrant
        4.0,                        # β2: 2nd beta dist param for kp entrant

        0.75, # From rer98          # cu: capacity utilization for cp

        0.02,                       # ϵ: minimum desired wage increase rate
        0.1,                        # max_g_wᴼ: max growth rate of offered wages
        0.5,                        # Kg_max: maximum capital growth rate

        0.5,                        # ωD: memory parameter cp demand estimation
        0.5,                        # ωQ: memory parameter cp quantity estimation
        0.5,                        # ωL: memory parameter cp labor supply estimation
        0.5,                        # ωW: memory parameter cp offered wage parameter

        0.95,                       # α_cp: parameter controlling MPC of consumers
        0.7,                        # c_L_max
        1000,                       # a_σ
        30,                         # b_σ

        0.05,                       # ψ_E: chance of employed worker looking for a better paying job
        0.75,                       # ψ_Q: chance of household switching away from cp when demand constrained
        0.75,                       # ψ_P: chance of household switching to cp with better price

        # 0.4                         # Nᵈ_share: share of expected demand cp wants to have in inventory
    )

    # Change parameters if needed before returning.
    if changed_params !== nothing

        for (key, new_param) in changed_params
            # println("$key, $new_param")
            setproperty!(global_param, Symbol(key), new_param)
        end
    end


    # df = DataFrame(Dict(x=>getfield(global_param, x) for x in fieldnames(GlobalParam)))
    # CSV.write("parameters/param_global_default.csv", df)

    return global_param
end