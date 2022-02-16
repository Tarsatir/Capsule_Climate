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