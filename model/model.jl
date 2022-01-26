
mutable struct All_Agents
    households :: Array{AbstractAgent}
    capital_good_producers :: Array{AbstractAgent}
    consumer_good_producers :: Array{AbstractAgent}
    capital_good_euclidian_matrix :: Array
end

struct GlobalParam
    # F1 :: Int                       # number of kp
    # F2 :: Int                       # number of cp
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
    # r :: Float64                    # interest rate
    φ1 :: Float64                   # 1st Uniform dist support, cp entrant cap
    φ2 :: Float64                   # 2nd Uniform dist support, cp entrant cap
    φ3 :: Float64                   # 1st Uniform dist support, cp entrant liq
    φ4 :: Float64                   # 2nd Uniform dist support, cp entrant liq
    α2 :: Float64                   # 1st beta dist param for kp entrant
    β2 :: Float64                   # 2nd beta dist param for kp entrant
    # ψ1 :: Float64
    # ψ2 :: Float64
    # ψ3 :: Float64                 # param changed in current model
    # tr :: Float64
    # ϕ :: Float64
    cu :: Float16

    ϵ :: Float64                    # minimum desired wage increase rate

    ωD :: Float64                   # memory parameter cp demand estimation
    ωQ :: Float64                   # memory parameter cp quantity estimation
    ωL :: Float64                   # memory parameter cp labor supply estimation
end

function initialize_global_params()
    global_param = GlobalParam(
        # 5, # TEMP
        # 20, # TEMP
        0.04,
        0.5,
        0.3,
        3.0,
        3.0,
        -0.15,
        0.15,
        0.5,
        0.04,
        0.1,
        3,
        20,
        0.04,
        1.0,
        1.0,
        2.0,
        # 0.01,
        0.1,
        0.9,
        0.1,
        0.9,
        2.0,
        4.0,
        # 1.0,
        # 0.0,
        # 0.0,
        # 0.1,
        # 0.4,
        0.75, # From rer98

        0.02,                       # ϵ: minimum desired wage increase rate

        0.5,                        # ωD: memory parameter cp demand estimation
        0.5,                        # ωQ: memory parameter cp quantity estimation
        0.5                         # ωL: memory parameter cp labor supply estimation
    )
    return global_param
end

function initialize_model()
    model_struct = Model([], [], [])
    return model_struct
end

function get_capgood_euclidian(all_agents, n_captlgood)
    """
    Generates a matrix with Euclidian distances of the 
    technological levels A and B. 
    """
    distance_matrix = zeros((n_captlgood, n_captlgood))

    all_kp = all_agents.capital_good_producers
    
    for i in 1:n_captlgood
        for j in i:n_captlgood

            # get current values for A and B of both producers
            A1 = all_kp[i].A[end]
            A2 = all_kp[j].A[end]
            B1 = all_kp[i].B[end]
            B2 = all_kp[j].B[end]
            
            distance = sqrt((A1-A2)^2 + (B1-B2)^2)
            if (i==j)
                distance = Inf
            end
            distance_matrix[i,j] = distance
            distance_matrix[j,i] = distance
        end
    end 
    all_agents.capital_good_euclidian_matrix = distance_matrix
end
