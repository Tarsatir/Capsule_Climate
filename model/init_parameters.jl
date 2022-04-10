@with_kw struct InitParam
    # Agent counts
    n_kp::Int = 20                        # number of kp
    n_cp::Int = 200                       # number of cp
    n_hh::Int = 2500                      # number of hh

    # Init params of kp
    n_init_emp_kp::Int = 4                # number of employees of kp at init
    A_LP_0::Float64 = 1.0                 # initial productivity level A_LP
    A_EE_0::Float64 = 1.0                 # initial productivity level A_EE
    A_EF_0::Float64 = 1.0                 # initial productivity level A_EF
    B_0::Float64 = 1.0                    # initial productivity level B

    # Init params of cp
    n_init_emp_cp::Int = 12               # number of employees of cp at init
    n_machines_init::Int = 50             # number of machines of cp at init

    # Init params of hh
    n_bp_hh::Int = 7                      # number of bp of hh (also min amount)
    n_lp_hh::Int = 7                      # number of lp of hh (also min amount)
end