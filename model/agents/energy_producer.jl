@with_kw mutable struct EnergyProducer
    T::Int=T                                    # Total number of iterations

    Dₑ::Vector{Float64} = zeros(Float64, T)     # Demand for energy units over time

    # Prices, cost and investments
    μₑ::Float64 = 0.2                           # Markup to determine price
    Πₑ::Vector{Float64} = zeros(Float64, T)     # Profits over time
    pₑ::Vector{Float64} = zeros(Float64, T)     # Price of energy over time
    PCₑ::Vector{Float64} = zeros(Float64, T)    # Cost of generating Dₑ(t) units of energy over time
    ICₑ::Vector{Float64} = zeros(Float64, T)    # Expansion and replacement investments over time
    RDₑ::Vector{Float64} = zeros(Float64, T)    # R&D expenditure over time

    # Owned power plants
    green_portfolio::Dict{Int, PowerPlant} = Dict()   # Portfolio of all green powerplants
    capacity_green::Float64 = 0.0               # capacity of set of green powerplants
    dirty_portfolio::SortedDict{Int, PowerPlant} = SortedDict() # Portfolio of all dirty powerplants
    capacity_dirty::Float64 = 0.0               # capacity of set of dirty powerplants
    infra_marg::Vector{Int} = []                # Set of ids of infra-marginal plants
end


"""
Production process of ep
"""
function produce_energy_ep!(
    ep::EnergyProducer,
    t::Int
    )


end


"""
Investment process of ep
"""
function plan_investments_ep!(
    ep::EnergyProducer, 
    t::Int
    )

    
end

"""
Computes the profits resulting from last periods production.
    Lamperti (2018), eq 10.
"""
function compute_Πₑ_ep!(
    ep::EnergyProducer, 
    t::Int
    )

    ep.Πₑ[t] = ep.pₑ[t] * ep.Dₑ[t] - ep.PCₑ[t] - ep.ICₑ[t] - ep.RDₑ[t]
end


"""
Computes cost of production PCₑ of infra marginal power plants
    Lamperti (2018), eq 12.
"""
function compute_PCₑ_ep!(
    ep::EnergyProducer,
    t::Int
    )

    dirty_cost = 0.0

    for pp_id in ep.infra_marg
        if haskey(ep.dirty_portfolio, pp_id)
            dirty_cost += ep.dirty_portfolio[pp].freq * ep.dirty_portfolio[pp].c * ep.dirty_portfolio[pp].Aᵗ
        end
    end

    ep.PCₑ[t] = dirty_cost
end


"""
Computes the price for each sold energy unit
"""
function compute_pₑ_ep!(
    ep::EnergyProducer, 
    t::Int
    )

    c̄ = ep.dirty_portfolio[ep.infra_marg[end]].c
    ep.pₑ[t] = ep.Dₑ[t] <= ep.capacity_green ? ep.μₑ : c̄ + ep.μₑ
end



"""
Lets ep choose power plants to produce energy demand with
"""
function choose_powerplants_ep!(
    ep::EnergyProducer,
    t::Int
    )

    # Check if all production can be done using green tech, if not, compute cost
    # of production using dirty tech
    if ep.Dₑ[t] < ep.capacity_green
        ep.infra_marg = collect(keys(ep.green_portfolio))
    else

        ep.infra_marg = collect(keys(ep.green_portfolio))
        total_capacity = ep.capacity_green

        # Sort dirty portfolio as to take low cost power plants first
        sort!(ep.dirty_portfolio, by = pp -> pp.c)

        for (i, powerplant) in ep.dirty_portfolio
            push!(ep.infra_marg, i)
            total_capacity += powerplant.capacity
            if total_capacity >= ep.Dₑ
                break
            end
        end
    end
end