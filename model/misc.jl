mutable struct Machine
    A :: Float64                # labor productivity machine
    freq :: Float64             # freq machine owned by cp
    age :: Float64              # age of machine
end

function initialize_machine(
    freq::Float64,
    η::Int=0, 
    A::Float64=1.0
    )

    machine_struct = Machine(
        A,                      # A: labor productivity machine
        freq,                   # freq: freq machine owned by cp
        sample(0:η)             # age: age of machine
    )
    return machine_struct
end


"""
Generates a matrix with Euclidian distances of the 
technological levels A and B. 
"""
function get_capgood_euclidian(
    all_kp::Vector{Int}, 
    model::ABM
    )::Array

    n_captlgood = length(all_kp)
    distance_matrix = zeros((n_captlgood, n_captlgood))

    for i in 1:n_captlgood
        for j in i:n_captlgood

            # get current values for A and B of both producers
            A1 = model[all_kp[i]].A[end]
            A2 = model[all_kp[j]].A[end]
            B1 = model[all_kp[i]].B[end]
            B2 = model[all_kp[j]].B[end]
            
            distance = sqrt((A1-A2)^2 + (B1-B2)^2)
            if (i==j)
                distance = Inf
            end
            distance_matrix[i,j] = distance
            distance_matrix[j,i] = distance
        end
    end 
    return distance_matrix
end