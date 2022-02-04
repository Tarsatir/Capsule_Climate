"""
    per_type(type, shuffle_agents::Bool)
A scheduler that returns multiple schedulers, each for a specific subset 
`shuffle_agents = true` randomizes the order of agents within each group.
"""
function per_type(shuffle_agents::Bool, model::ABM)

    all_hh = Vector{Int}()
    all_cp = Vector{Int}()
    all_kp = Vector{Int}()
    all_bp = Vector{Int}()
    all_lp = Vector{Int}()
    all_p = Vector{Int}()

    for agent in allagents(model)
        if typeof(agent) == Household
            push!(all_hh, agent.id)
        elseif typeof(agent) == ConsumerGoodProducer
            push!(all_cp, agent.id)
            push!(all_p, agent.id)
            if agent.type_good == "Basic"
                push!(all_bp, agent.id)
            else
                push!(all_lp, agent.id)
            end
        elseif typeof(agent) == CapitalGoodProducer
            push!(all_kp, agent.id)
            push!(all_p, agent.id)
        end
    end

    if shuffle_agents
        for set in [all_hh, all_cp, all_kp, all_bp, all_lp]
            shuffle!(set)
        end
    end

    return all_hh, all_cp, all_kp, all_bp, all_lp, all_p

end