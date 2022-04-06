struct InitParam
    # Agent counts
    n_kp :: Int                             # number of kp
    n_cp :: Int                             # number of cp
    n_hh :: Int                             # number of hh

    # Init params of kp
    n_init_emp_kp :: Int                    # number of employees of kp at init
    A_0 :: Float64                          # initial productivity level A
    B_0 :: Float64                          # initial productivity level B

    # Init params of cp
    n_init_emp_cp :: Int                    # number of employees of cp at init
    n_machines_init :: Int                  # number of machines of cp at init

    # Init params of hh
    n_bp_hh :: Int                          # number of bp of hh (also min amount)
    n_lp_hh :: Int                          # number of lp of hh (also min amount)
end

function initialize_init_params()
    init_params = InitParam(
        # Agent counts
        20,                                 # n_kp: number of kp
        200,                                # n_cp: number of cp
        2500,                               # n_hh: number of hh

        # Init params of kp
        4,                                  # n_init_emp_kp: number of employees of kp at init
        2.5,                                # A_0: initial productivity level A
        2.5,                                # B_0: initial productivity level B

        # Init params of cp
        12,                                 # n_init_emp_cp: number of employees of cp at init
        50,                                 # n_machines_init: number of machines of cp at init

        # Init params of hh
        7,                                  # n_bp_hh: number of bp of hh (also min amount)
        7                                   # n_lp_hh: number of lp of hh (also min amount)
    )
    return init_params
end