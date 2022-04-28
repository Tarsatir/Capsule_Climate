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

    weights::Matrix{Float64} = rand(n_hh, n_cp)
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

    for (i,hh_id) ∈ enumerate(all_hh)
        cm_dat.all_C[i] = model[hh_id].C
        for (j,cp_id) ∈ enumerate(all_cp)
            cm_dat.weights[i,j] = 1 / model[cp_id].p[end]^2
            cm_dat.all_N[j] = model[cp_id].N_goods * model[cp_id].p[end]
        end
    end

    # for i in 1:cm_dat.n_hh

    #     cm_dat.weights_sum[i] = 0.0
    #     cm_dat.sold_per_hh[i] = 0.0
    #     cm_dat.sold_per_hh_round[i] = 0.0

    #     for j in 1:cm_dat.n_cp
    #         cm_dat.true_D[i,j] = 0.0
    #         cm_dat.transactions[i,j] = 0.0

    #         cm_dat.demand_per_cp[j] = 0.0
    #         cm_dat.sold_per_cp[j] = 0.0
    #         cm_dat.sold_per_cp_round[j] = 0.0
    #     end
    # end

    cm_dat.true_D .= 0.0
    cm_dat.weights_sum .= 0.0
    cm_dat.transactions .= 0.0
    cm_dat.demand_per_cp .= 0.0

    cm_dat.sold_per_hh .= 0.0
    cm_dat.sold_per_hh_round .= 0.0
    cm_dat.sold_per_cp .= 0.0
    cm_dat.sold_per_cp_round .= 0.0
end