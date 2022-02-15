mutable struct Balance
    # assets
    N :: Float64            # inventories
    K :: Float64            # capital
    NW :: Float64           # liquid assets

    # liabilities
    Deb :: Float64          # debt
    EQ :: Float64           # equity
end

mutable struct GovCurrentAccount
    # Revenues
    Rev_τᴵ :: Vector{Float64}            # hist revenues of income tax
    Rev_τˢ :: Vector{Float64}            # hist revenues of sales tax
    Rev_τᴾ :: Vector{Float64}            # hist revenues of profit tax
    Rev_τᴱ :: Vector{Float64}            # hist revenues of energy tax
    Rev_τᶜ :: Vector{Float64}            # hist revenues of emission tax

    # Expenditures
    Exp_UB :: Vector{Float64}            # hist expenditures on unemployment benefits
    Exp_Sub :: Vector{Float64}           # hist expenditures on subsidies

end

# function initialize_govcurrentaccount()
#     gov_curr_acc = GovCurrentAccount(
#         [],
#         [],
#         [],
#         [],
#         [],
#         [],
#         []
#     )
#     return gov_curr_acc
# end

mutable struct Machine
    A :: Float64                # labor productivity machine
    c :: Float64                # cost to produce machine
    freq :: Float64             # freq machine owned by cp
    age :: Float64              # age of machine
end

function initialize_machine(η=0)
    machine_struct = Machine(
        1,                      # A: labor productivity machine
        0,                      # c: cost to produce machine
        40,                     # freq: freq machine owned by cp
        rand(0:η)               # age: age of machine
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