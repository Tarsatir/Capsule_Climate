@with_kw struct InitParam
    
    # Agent counts
    n_kp::Int64 = 20                      # number of kp
    n_cp::Int64 = 200                     # number of cp
    n_hh::Int64 = 2500                    # number of hh

    # Init rates
    init_unempl_rate::Float64 = 0.05      # initial unemployment rate

    # Init params of kp
    # n_init_emp_kp::Int64 = 10             # number of employees of new kp at init
    A_LP_0::Float64 = 1.0                 # initial productivity level A_LP
    A_EE_0::Float64 = 1.0                 # initial productivity level A_EE
    A_EF_0::Float64 = 1.0                 # initial productivity level A_EF
    B_LP_0::Float64 = 1.0                 # initial productivity level B_LP
    B_EE_0::Float64 = 1.0                 # initial productivity level B_EE
    B_EF_0::Float64 = 1.0                 # initial productivity level B_EF

    # Init params of cp
    # n_init_emp_cp::Int = 10               # number of employees of cp at init
    n_machines_init::Int64 = 40           # number of machines of cp at init

    # Init params of hh
    n_cp_hh::Int64 = 25                   # number of cp of hh (also min amount)
    skill_mean::Float64 = 0.              # mean value for lognormal skill distribution
    skill_var::Float64 = 0.75             # variance value for lognormal skill distribution
    βmin::Float64 = 0.7                   # minimum value discount factor
    βmax::Float64 = 1.0                   # maximum value discount factor

    # Init params of ep
    n_powerplants_init::Int64 = 300_000   # number of unit of power plants in ep
    frac_green::Float64 = 0.1             # fraction of initial power plants that are green
    markup_ep::Float64 = 0.01             # markup rate energy producer
    Aᵀ_0::Float64 = 1.0                   # initial thermal efficiency
    emᵀ_0::Float64 = 1.0                  # initial emission level
    IC_g_0::Float64 = 12                  # initial fixed costs of green plant investments
end