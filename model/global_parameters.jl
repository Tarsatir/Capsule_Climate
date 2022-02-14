struct GlobalParam
    ν :: Float64                    # R&D inv propensity
    ξ :: Float64                    # R&D allocation to IN
    ζ :: Float64                    # firm search capabilities
    α1 :: Float64                   # 1st beta dist param for IN
    β1 :: Float64                   # 2nd beta dist param for IN
    κ_lower :: Float64              # 1st beta dist support
    κ_upper :: Float64              # 2nd beta dist support
    γ :: Float64                    # new custommer sample parameter
    μ1 :: Float64                   # kp markup rule
    ι :: Float64                    # desired inventories
    b :: Int                        # payback period
    η :: Int                        # physical scrapping age
    υ :: Float64                    # markup coefficient
    ω :: Float64                    # competetiveness weights
    χ :: Float64                    # replicator dynamics coeff
    Λ :: Float64                    # max debt/sales ratio
    φ1 :: Float64                   # 1st Uniform dist support, cp entrant cap
    φ2 :: Float64                   # 2nd Uniform dist support, cp entrant cap
    φ3 :: Float64                   # 1st Uniform dist support, cp entrant liq
    φ4 :: Float64                   # 2nd Uniform dist support, cp entrant liq
    α2 :: Float64                   # 1st beta dist param for kp entrant
    β2 :: Float64                   # 2nd beta dist param for kp entrant
    cu :: Float64                   # capacity utilization for cp
    ϵ :: Float64                    # minimum desired wage increase rate
    ωD :: Float64                   # memory parameter cp demand estimation
    ωQ :: Float64                   # memory parameter cp quantity estimation
    ωL :: Float64                   # memory parameter cp labor supply estimation
    α_cp :: Float64                 # parameter controlling MPC of consumers
    c_L_max :: Float64              # maximum share consumed on luxury goods
    a_σ :: Float64                  # 1st parameter governing logistic function
    b_σ :: Float64                  # 2nd parameter governing logistic function
end


function initialize_global_params()
    global_param = GlobalParam(
        0.04,                       # ν: R&D inv propensity
        0.5,                        # ξ: R&D allocation to IN
        0.3,                        # ζ: firm search capabilities
        3.0,                        # α1: 1st beta dist param for IN
        3.0,                        # β1: 2nd beta dist param for IN
        -0.15,                      # κ_lower: 1st beta dist support
        0.15,                       # κ_upper: 2nd beta dist support
        0.5,                        # γ: new custommer sample parameter
        0.2,                        # μ1: kp markup rule
        0.1,                        # ι: desired inventories
        3,                          # b: payback period
        20,                         # η: physical scrapping age
        0.04,                       # υ: markup coefficient
        1.0,                        # ω: competetiveness weights
        1.0,                        # χ: replicator dynamics coeff
        2.0,                        # Λ: max debt/sales ratio
        0.1,                        # φ1: 1st Uniform dist support, cp entrant cap
        0.9,                        # φ2: 2nd Uniform dist support, cp entrant cap
        0.1,                        # φ3: 1st Uniform dist support, cp entrant liq
        0.9,                        # φ4: 2nd Uniform dist support, cp entrant liq
        2.0,                        # α2: 1st beta dist param for kp entrant
        4.0,                        # β2: 2nd beta dist param for kp entrant
        0.75, # From rer98          # cu: capacity utilization for cp
        0.02,                       # ϵ: minimum desired wage increase rate
        0.1,                        # ωD: memory parameter cp demand estimation
        0.1,                        # ωQ: memory parameter cp quantity estimation
        0.1,                        # ωL: memory parameter cp labor supply estimation
        0.9,                        # α_cp: parameter controlling MPC of consumers
        0.7,                        # c_L_max
        100,                        # a_σ
        3                           # b_σ
    )
    return global_param
end