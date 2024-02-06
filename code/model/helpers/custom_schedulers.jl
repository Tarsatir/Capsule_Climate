"""
    schedule_per_type(type, shuffle_agents::Bool)
A scheduler that returns multiple schedulers, each for a specific subset 
`shuffle_agents = true` randomizes the order of agents within each group.
"""
function schedule_per_type(
    model::ABM;
    shuffle_agents=true::Bool
    )

    if shuffle_agents
        shuffle!(model.all_hh)
        shuffle!(model.all_cp)
        shuffle!(model.all_kp)
        shuffle!(model.all_p)
    end

    # all_hh = Vector{Int}()
    # all_cp = Vector{Int}()
    # all_kp = Vector{Int}()
    # all_p = Vector{Int}()

    # for agent in allagents(model)
    #     if typeof(agent) == Household
    #         push!(all_hh, agent.id)
    #     elseif typeof(agent) == ConsumerGoodProducer
    #         push!(all_cp, agent.id)
    #         push!(all_p, agent.id)
    #     elseif typeof(agent) == CapitalGoodProducer
    #         push!(all_kp, agent.id)
    #         push!(all_p, agent.id)
    #     end
    # end

    # if shuffle_agents
    #     for set in [all_hh, all_cp, all_kp, all_p]
    #         shuffle!(set)
    #     end
    # end

    return model.all_hh, model.all_cp, model.all_kp, model.all_p
end