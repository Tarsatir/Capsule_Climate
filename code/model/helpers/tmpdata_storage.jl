

"""
CMDATA - DATA STRUCT FOR CONSUMER MARKET DATA
"""

"""
Mutable struct that holds the data structures used to save data from the consumer
    market process.
"""
@with_kw mutable struct CMData
    n_hh::Int
    n_cp::Int

    true_D::Matrix{Float64} = spzeros(Float64, n_hh, n_cp)
    all_C::Vector{Float64} = zeros(Float64, n_hh)
    all_N::Vector{Float64} = zeros(Float64, n_cp)

    sold_per_hh::Vector{Float64} = spzeros(Float64, n_hh)
    sold_per_hh_round::Vector{Float64} = zeros(Float64, n_hh)
    sold_per_cp::Vector{Float64} = zeros(Float64, n_cp)
    sold_per_cp_round::Vector{Float64} = zeros(Float64, n_cp)

    weights::Matrix{Float64} = spzeros(Float64, n_hh, n_cp)
    weights_sum::Vector{Float64} = zeros(Float64, n_hh)

    transactions::Matrix{Float64} = spzeros(Float64, n_hh, n_cp)
    frac_sellable::Vector{Float64} = ones(Float64, n_cp)
    C_spread::Matrix{Float64} = spzeros(Float64, n_hh, n_cp)
    demand_per_cp::Vector{Float64} = zeros(Float64, n_cp)
end


"""
Resets fields in cm data struct before cm market process is initiated
"""
function reset_matrices_cp!(
    cmdata::CMData,
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    model::ABM
    )

    # Set to order of small to large id (minimum(all_cp):max(all_cp)
    all_p = map(cp_id -> model[cp_id].p[end], minimum(all_cp):maximum(all_cp))
    all_N_goods = map(cp_id -> model[cp_id].N_goods, minimum(all_cp):maximum(all_cp))

    @inbounds for (i,hh_id) ∈ enumerate(minimum(all_hh):maximum(all_hh))

        cmdata.all_C[i] = model[hh_id].C
        cmdata.weights[i,:] .= 0.0

        for cp_id in model[hh_id].cp

            # cp are initiated after hh and cp_id will thus correspond to col + len(hh), 
            # so subtract len(hh) to get index
            j = cp_id - length(all_hh)

            cmdata.all_N[j] = all_N_goods[j] * all_p[j]
            cmdata.weights[i,j] = all_p[j] ^ -1
        end
    end

    cmdata.transactions .= 0.0
end


"""
GINIDATA - DATA STRUCT FOR CONSUMER MARKET DATA
"""

"""
Struct that holds intermediate data for GINI computations
"""
@with_kw mutable struct GINIData{D}
    I::Matrix{Float64} = Matrix{Float64}(undef, D, D)
    W::Matrix{Float64} = Matrix{Float64}(undef, D, D)
end