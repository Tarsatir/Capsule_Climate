mutable struct ConsumerMarket
    demanding_hh :: Array{AbstractAgent}
    supplying_bp :: Array{AbstractAgent}
    supplying_lg :: Array{AbstractAgent}
end


function initialize_consumermarket()
    consumermarket_struct = ConsumerMarket(
        [],
        [],
        []
    )
    return consumermarket_struct
end


function consumermarket_process!(consumermarket_struct, all_hh :: Array, all_bp :: Array, all_lp :: Array, gov_struct)

    # TODO: put this as a parameter somewhere
    n_rounds = 1

    demanding_hh = copy(all_hh)
    supplying_bp = copy(all_bp)
    supplying_lp = copy(all_lp)

    # add zero to all hist demand for this time step
    for p in vcat(all_bp, all_lp)
        push!(p.D, 0)
    end

    # loop over all demanding households for n rounds and match producers
    for round in 1:n_rounds

        # TODO: find a way so cons package not determined all over again, otherwise 
        # consumers buy too much bg or lg
        
        for hh in demanding_hh

            if length(supplying_lp) == 0 || length(supplying_bp) == 0
                return
            end
            
            # pick bp and lp 
            pick_cp_hh!(hh, supplying_bp, supplying_lp)

            # set budget based on bp and lp, decide consumption amount
            N_B, N_L = set_cons_package_hh!(hh, gov_struct.τˢ)

            # transact goods if enough available
            bg_satisfied =  transact_cp!(hh.bg, hh, N_B)   
            lg_satisfied =  transact_cp!(hh.lg, hh, N_L)

            # implies both bg and lg have inventory left
            if bg_satisfied && lg_satisfied
                filter!(h -> h ≠ hh, demanding_hh)
            elseif bg_satisfied && ~(lg_satisfied)
                filter!(lp -> lp ≠ hh.lg, supplying_lp)
            elseif lg_satisfied && ~(bg_satisfied)
                filter!(bp -> bp ≠ hh.bg, supplying_bp)
            elseif ~(bg_satisfied && lg_satisfied)
                filter!(bp -> bp ≠ hh.bg, supplying_bp)
                filter!(lp -> lp ≠ hh.lg, supplying_lp)
            end
        end
    end

    # match all remaining households and producers randomly



end