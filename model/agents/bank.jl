mutable struct Bank 
    Deposits::Dict{}               # deposits at bank
    Debts::Dict{}                  # debts at bank
end

function initialize_bank(
    all_p::Vector{Int},
    model::ABM
    )::Bank

    # Add deposits for all agents
    all_deposits = Dict()
    for agent in allagents(model)
        all_deposits[agent.id] = 0.0
    end

    # Add debts (only for producers)
    all_debts = Dict()
    for p_id in all_p
        all_debts[p_id] = 0.0
    end

    # Initialize bank struct
    bank_struct = Bank(
        all_deposits,
        all_debts
    )
    return bank_struct

end

