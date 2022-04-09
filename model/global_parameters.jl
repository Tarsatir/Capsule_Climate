@with_kw mutable struct GlobalParam
    # Determine technical innovation process
    ν::Float64 = 0.04               # R&D inv propensity
    ξ::Float64 = 0.5                # R&D allocation to IN
    ζ::Float64 = 0.3                # firm search capabilities
    α1::Float64 = 3.0               # 1st beta dist param for IN
    β1::Float64 = 3.0               # 2nd beta dist param for IN
    κ_lower::Float64 = -0.02        # 1st beta dist support
    κ_upper::Float64 = 0.02         # 2nd beta dist support

    γ::Float64 = 0.5                # new custommer sample parameter
    μ1::Float64 = 0.2               # kp markup rule
    r::Float64 = 0.0                # interest rate
    ι::Float64 = 0.1                # desired inventories
    b::Int = 3                      # payback period
    η::Int = 20                     # physical scrapping age
    Λ::Float64 = 2.0                # max debt/sales ratio, regular debt

    # Determine entrant composition
    φ1::Float64 = 0.1               # 1st Uniform dist support, cp entrant cap
    φ2::Float64 = 0.9               # 2nd Uniform dist support, cp entrant cap
    φ3::Float64 = 0.1               # 1st Uniform dist support, cp entrant liq
    φ4::Float64 = 0.9               # 2nd Uniform dist support, cp entrant liq
    α2::Float64 = 2.0               # 1st beta dist param for kp entrant
    β2::Float64 = 4.0               # 2nd beta dist param for kp entrant
    φ5::Float64 = -0.04             # 1st Beta dist support for kp entrant tech
    φ6::Float64 = 0.02              # 2nd Beta dist support for kp entrant tech

    cu::Float64 = 0.75              # capacity utilization for cp
    max_NW_ratio::Float64 = 2.0     # maximum ratio p can have monthly expenses in NW
    ϵ::Float64 = 0.02               # minimum desired wage increase rate
    max_g_wᴼ::Float64 = 0.1         # max growth rate of offered wages
    Kg_max::Float64 = 0.5           # maximum capital growth rate

    # Determine expectation updating cp
    ω::Float64 = 0.5                # memory parameter adaptive updating rules

    # Determine household consumption
    α_cp::Float64 = 0.9             # parameter controlling MPC of consumers
    c_L_max::Float64 = 0.7          # maximum share consumed on luxury goods
    a_σ::Float64 = 1000             # 1st parameter governing logistic function
    b_σ::Float64 = 30               # 2nd parameter governing logistic function

    # Determine household switching
    ψ_E::Float64 = 0.15             # chance of employed worker looking for a better paying job
    ψ_Q::Float64 = 0.75             # chance of household switching away from cp when demand constrained
    ψ_P::Float64 = 0.75             # chance of household switching to cp with better price

    freq_per_machine::Int = 50      # capital units per machine

    n_cons_market_days::Int = 4     # number of days in the consumer market process

    t_wait::Int = 4                 # number of time periods producers are not allowed to go bankrupt

    fordist_lm::Bool = false        # determines whether the labormarket is Fordist or competetive
end


function initialize_global_params(
    labormarket_is_fordist::Bool,
    changed_params
    )

    global_param = GlobalParam()
    # global_param = GlobalParam(
    #     0.04,                       # ν: R&D inv propensity
    #     0.5,                        # ξ: R&D allocation to IN
    #     0.3,                        # ζ: firm search capabilities
    #     3.0,                        # α1: 1st beta dist param for IN
    #     3.0,                        # β1: 2nd beta dist param for IN
    #     -0.02,                      # κ_lower: 1st beta dist support
    #     0.02,                       # κ_upper: 2nd beta dist support

    #     0.5,                        # γ: new custommer sample parameter
    #     0.2,                        # μ1: kp markup rule
    #     0.0,                        # r: (annual) interest rate
    #     0.1,                        # ι: desired inventories
    #     3,                          # b: payback period
    #     20,                         # η: physical scrapping age
    #     2.0,                        # Λ: max debt/sales ratio regular debt

    #     0.1,                        # φ1: 1st Uniform dist support, cp entrant cap
    #     0.9,                        # φ2: 2nd Uniform dist support, cp entrant cap
    #     0.1,                        # φ3: 1st Uniform dist support, cp entrant liq
    #     0.9,                        # φ4: 2nd Uniform dist support, cp entrant liq
    #     2.0,                        # α2: 1st beta dist param for kp entrant
    #     4.0,                        # β2: 2nd beta dist param for kp entrant

    #     0.75, # From rer98          # cu: capacity utilization for cp

    #     0.02,                       # ϵ: minimum desired wage increase rate
    #     0.1,                        # max_g_wᴼ: max growth rate of offered wages
    #     0.5,                        # Kg_max: maximum capital growth rate

    #     0.5,                        # ω: memory parameter adaptive updating rules

    #     0.9,                       # α_cp: parameter controlling MPC of consumers
    #     0.7,                        # c_L_max
    #     1000,                       # a_σ
    #     30,                         # b_σ

    #     0.15,                       # ψ_E: chance of employed worker looking for a better paying job
    #     0.75,                       # ψ_Q: chance of household switching away from cp when demand constrained
    #     0.75,                       # ψ_P: chance of household switching to cp with better price

    #     50,                         # capital units per machine

    #     4,                           # number of days in consumer market proces
    #     labormarket_is_fordist
    # )

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