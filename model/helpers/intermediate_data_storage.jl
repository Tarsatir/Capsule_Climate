"""
Mutable struct that holds the data structures used to save data from the consumer
    market process.
"""
@with_kw mutable struct CMData
    n_hh::Int
    n_cp::Int

    true_D::Matrix{Float64} = zeros(Float64, n_hh, n_cp)
    all_C::Vector{Float64} = zeros(Float64, n_hh)
    all_N::Vector{Float64} = zeros(Float64, n_cp)

    sold_per_hh::Vector{Float64} = zeros(Float64, n_hh)
    sold_per_hh_round::Vector{Float64} = zeros(Float64, n_hh)
    sold_per_cp::Vector{Float64} = zeros(Float64, n_cp)
    sold_per_cp_round::Vector{Float64} = zeros(Float64, n_cp)

    weights::Matrix{Float64} = zeros(Float64, n_hh, n_cp)
    weights_sum::Vector{Float64} = zeros(Float64, n_hh)

    transactions::Matrix{Float64} = zeros(Float64, n_hh, n_cp)
    frac_sellable::Vector{Float64} = ones(Float64, n_cp)
    C_spread::Matrix{Float64} = zeros(Float64, n_hh, n_cp)
    demand_per_cp::Vector{Float64} = zeros(Float64, n_cp)
end


"""
Resets fields in cm data struct before cm market process is initiated
"""
function reset_matrices_cp!(
    cm_dat::CMData,
    all_hh::Vector{Int},
    all_cp::Vector{Int},
    model::ABM
    )

    # Set to order of small to large id (minimum(all_cp):max(all_cp)
    all_p = map(cp_id -> model[cp_id].p[end], minimum(all_cp):maximum(all_cp))
    all_N_goods = map(cp_id -> model[cp_id].N_goods, minimum(all_cp):max(all_cp))

    @inbounds for (i,hh_id) ∈ enumerate(minimum(all_hh):max(all_hh))

        cm_dat.all_C[i] = model[hh_id].C
        cm_dat.weights[i,:] .= 0.0

        # for j in indexin(model[hh_id].cp, all_cp)
        for cp_id in model[hh_id].cp

            # cp are initiated after hh and cp_id will thus correspond to col + len(hh), 
            # so subtract len(hh) to get index
            j = cp_id - length(all_hh)

        # @inbounds for (j,cp_id) ∈ enumerate(all_cp)

            cm_dat.all_N[j] = all_N_goods[j] * all_p[j]
            cm_dat.weights[i,j] = 1 / all_p[j]^2

            # if cp_id ∈ model[hh_id].cp
                # cm_dat.weights[i,j] = 1 / all_p[j]^2
            # else
                # cm_dat.weights[i,j] = 0.0
            # end

        end
    end

    # cm_dat.true_D .= 0.0
    # cm_dat.weights_sum .= 0.0
    cm_dat.transactions .= 0.0
    # cm_dat.demand_per_cp .= 0.0

    # cm_dat.sold_per_hh .= 0.0
    # cm_dat.sold_per_hh_round .= 0.0
    # cm_dat.sold_per_cp .= 0.0
    # cm_dat.sold_per_cp_round .= 0.0
end