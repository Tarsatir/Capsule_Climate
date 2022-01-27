mutable struct ConsumerMarket
    demanding_hh :: Array{AbstractAgent}
    supplying_bg :: Array{AbstractAgent}
    supplying_lg :: Array{AbstractAgent}
end


function initialize_consumermarket()
    consumermarket_struct = ConsumerMarket(
        [],
        [],
        []
    )
end


function consumermarket_process!(all_hh, all_bp, all_lp, gov_struct)

    for hh in all_agents.all_hh
        pick_cp_hh(hh, all_agents.all_bp, all_agents.all_lp)
    end

end