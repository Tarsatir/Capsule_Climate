@with_kw struct InitParam
    # Agent counts
    n_kp::Int = 10                        # number of kp
    n_cp::Int = 200                       # number of cp
    n_hh::Int = 2500                      # number of hh

    # Init params of kp
    n_init_emp_kp::Int = 10                # number of employees of kp at init
    A_LP_0::Float64 = 1.0                 # initial productivity level A_LP
    A_EE_0::Float64 = 1.0                 # initial productivity level A_EE
    A_EF_0::Float64 = 1.0                 # initial productivity level A_EF
    B_LP_0::Float64 = 1.0                 # initial productivity level B_LP
    B_EE_0::Float64 = 1.0                 # initial productivity level B_EE
    B_EF_0::Float64 = 1.0                 # initial productivity level B_EF

    # Init params of cp
    n_init_emp_cp::Int = 10               # number of employees of cp at init
    n_machines_init::Int = 40             # number of machines of cp at init

    # Init params of hh
    n_cp_hh::Int = 25                     # number of cp of hh (also min amount)
    σ_hh_I::Float64 = 0.406               # sigma of lognormal distribution of income
    scale_hh_I::Float64 = 29191.575       # scale of lognormal distribution of income

    # Init params of ep
    n_powerplants_init::Int = 300_000     # number of unit of power plants in ep
    frac_green::Float64 = 0.1             # fraction of initial power plants that are green
    μₑ::Float64 = 0.01                     # Markup rate energy producer
    Aᵀ_0::Float64 = 1.0                   # initial thermal efficiency
    emᵀ_0::Float64 = 1.0                  # initial emission level
    IC_g_0::Float64 = 12                  # initial fixed costs of green plant investments
end