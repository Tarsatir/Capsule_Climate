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
            A1_LP = model[all_kp[i]].A_LP
            A2_LP = model[all_kp[j]].A_LP

            A1_EE = model[all_kp[i]].A_EE
            A2_EE = model[all_kp[j]].A_EE

            A1_EF = model[all_kp[i]].A_EF
            A2_EF = model[all_kp[j]].A_EF

            B1_LP = model[all_kp[i]].B_LP
            B2_LP = model[all_kp[j]].B_LP

            B1_EE = model[all_kp[i]].B_EE
            B2_EE = model[all_kp[j]].B_EE

            B1_EF = model[all_kp[i]].B_EF
            B2_EF = model[all_kp[j]].B_EF
            
            distance = sqrt((A1_LP-A2_LP)^2 + (A1_EE-A2_EE)^2 + (A1_EF-A2_EF)^2 + 
                            (B1_LP-B2_LP)^2 + (B1_EE-B2_EE)^2 + (B1_EF-B2_EF)^2)
            if (i==j)
                distance = Inf
            end
            distance_matrix[i,j] = distance
            distance_matrix[j,i] = distance
        end
    end 
    return distance_matrix
end