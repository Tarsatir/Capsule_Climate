# mutable struct ConsumerMarket
#     demanding_hh :: Vector{Int}
#     supplying_bp :: Vector{Int}
#     supplying_lg :: Vector{Int}
# end


# function initialize_consumermarket()
#     consumermarket_struct = ConsumerMarket(
#         Vector{Int}(),
#         Vector{Int}(),
#         Vector{Int}()
#     )
#     return consumermarket_struct
# end


function consumermarket_process!(
    # consumermarket_struct, 
    all_hh::Vector{Int},
    all_cp::Vector{Int}, 
    all_bp::Vector{Int}, 
    all_lp::Vector{Int}, 
    gov_struct,
    model::ABM
    )

    # TODO: put this as a parameter somewhere
    n_rounds = 1

    demanding_hh = copy(all_hh)
    supplying_bp = copy(all_bp)
    supplying_lp = copy(all_lp)

    # add zero to all hist demand for this time step
    for cp_id in all_cp
        push!(model[cp_id].D, 0)
    end

    # loop over all demanding households for n rounds and match producers
    for round in 1:n_rounds

        # TODO: find a way so cons package not determined all over again, otherwise 
        # consumers buy too much bg or lg
        
        for hh_id in demanding_hh

            hh = model[hh_id]

            if length(supplying_lp) == 0 || length(supplying_bp) == 0
                return
            end
            
            # pick bp and lp 
            pick_cp_hh!(hh, supplying_bp, supplying_lp)

            # set budget based on bp and lp, decide consumption amount
            N_B, N_L = set_cons_package_hh!(hh, gov_struct.τˢ, model)

            # println(N_B,  " ", N_L)

            # transact goods if enough available
            bg_satisfied = transact_cp!(model[hh.pref_bp_id], hh, N_B)   
            lg_satisfied = transact_cp!(model[hh.pref_lp_id], hh, N_L)

            # implies both bg and lg have inventory left
            if bg_satisfied && lg_satisfied
                filter!(h -> h ≠ hh, demanding_hh)
            elseif bg_satisfied && ~(lg_satisfied)
                filter!(lp -> lp ≠ hh.pref_lp_id, supplying_lp)
            elseif lg_satisfied && ~(bg_satisfied)
                filter!(bp -> bp ≠ hh.pref_bp_id, supplying_bp)
            elseif ~(bg_satisfied && lg_satisfied)
                filter!(bp -> bp ≠ hh.pref_bp_id, supplying_bp)
                filter!(lp -> lp ≠ hh.pref_lp_id, supplying_lp)
            end
        end
    end

    # println("Yeet ", length(all_bp[1].D))

    # match all remaining households and producers randomly


end