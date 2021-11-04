
mutable struct Model
    capital_good_producers :: Array
    consumer_good_producers :: Array
    capital_good_euclidian_matrix :: Array
end

struct GlobalParam
    F1 :: Int
    F2 :: Int
    ν :: Float64
    ξ :: Float64
    ζ :: Float64
    α1 :: Float64
    β1 :: Float64
    κ_lower :: Float64
    κ_upper :: Float64
    γ :: Float64
    μ1 :: Float64
    ι :: Float64
    b :: Int
    η :: Int
    υ :: Float64
    ω :: Float64
    χ :: Float64
    Λ :: Float64
    r :: Float64
    φ1 :: Float64
    φ2 :: Float64
    φ3 :: Float64
    φ4 :: Float64
    α2 :: Float64
    β2 :: Float64
    ψ1 :: Float64
    ψ2 :: Float64
    ψ3 :: Float64
    tr :: Float64
    ϕ :: Float64
    cu :: Float16
end

function initialize_global_params()
    global_param = GlobalParam(
        # 50,
        # 200,
        5, # TEMP
        20, # TEMP
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
        0.01,
        0.1,
        0.9,
        0.1,
        0.9,
        2.0,
        4.0,
        1.0,
        0.0,
        0.0,
        0.1,
        0.4
        0.75 # From rer98
    )
    return global_param
end

function initialize_model()
    model_struct = Model([], [], [])
    return model_struct
end

function get_capgood_euclidian(model_struct, global_param)
    """
    Generates a matrix with Euclidian distances of the 
    technological levels A and B. 
    """
    distance_matrix = zeros((global_param.F1, global_param.F1))
    
    for i in 1:global_param.F1
        for j in i:global_param.F1

            # get current values for A and B of both producers
            A1 = model_struct.capital_good_producers[i].A[end]
            A2 = model_struct.capital_good_producers[j].A[end]
            B1 = model_struct.capital_good_producers[i].B[end]
            B2 = model_struct.capital_good_producers[j].B[end]

            distance = sqrt((A1-A2)^2 + (B1-B2)^2)
            distance_matrix[i,j] = distance
            distance_matrix[j,i] = distance
        end
    end 
    # display(distance_matrix)
    model_struct.capital_good_euclidian_matrix = distance_matrix

end
