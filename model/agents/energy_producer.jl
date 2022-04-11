@with_kw mutable struct EnergyProducer
    T::Int=T                                    # Total number of iterations

    Dₑ::Vector{Float64} = zeros(Float64, T)     # Demand for energy units over time
    Q̄ₑ::Vector{Float64} = zeros(Float64, T)     # Maximum production of units over time

    # Prices, cost and investments
    μₑ::Float64 = 0.2                           # Markup to determine price
    Πₑ::Vector{Float64} = zeros(Float64, T)     # Profits over time
    pₑ::Vector{Float64} = zeros(Float64, T)     # Price of energy over time
    PCₑ::Vector{Float64} = zeros(Float64, T)    # Cost of generating Dₑ(t) units of energy over time
    ICₑ::Vector{Float64} = zeros(Float64, T)    # Expansion and replacement investments over time
    RDₑ::Vector{Float64} = zeros(Float64, T)    # R&D expenditure over time
    IN_g::Vector{Float64} = zeros(Float64, T)   # R&D spending allocated to innovation in green tech
    IN_d::Vector{Float64} = zeros(Float64, T)   # R&D spending allocated to innovation in dirty tech
    EIᵈ::Vector{Float64} = zeros(Float64, T)    # Desired type and amount of units of expansionary investments
    ECₑ::Vector{Float64} = zeros(Float64, T)    # Cost of spansionary investment

    # Technological parameters
    IC_g::Vector{Float64} = zeros(Float64, T)   # Fixed investment cost of the cheapest new green power plant
    Aᵀ_d::Vector{Float64} = zeros(Float64, T)   # Thermal efficiency of new power plants
    emᵀ_d::Vector{Float64} = zeros(Float64, T)  # Emissions of new power plants
    c_d::Vector{Float64} = zeros(Float64, T)    # Discounted production cost of the cheapest dirty plant 

    # Owned power plants
    green_portfolio::Vector{PowerPlant}         # Portfolio of all green powerplants
    green_capacity::Float64 = 0.0               # capacity of set of green powerplants
    dirty_portfolio::Vector{PowerPlant}         # Portfolio of all dirty powerplants
    dirty_capacity::Float64 = 0.0               # Capacity of set of dirty powerplants
    infra_marg::Vector{PowerPlant} = []         # Infra-marginal plants
    pp_tb_repl::Vector{PowerPlant} = []         # Power plants to be replaced
end


function initialize_energy_producer(
    T::Int,
    init_param::InitParam,
    global_param::GlobalParam
    )::EnergyProducer

    # Initialize power plants
    n_pp::Int = init_param.n_powerplants_init / global_param.freq_per_powerplant
    n_pp_green::Int = init_param.frac_green * n_pp

    # Initialize green power plants
    green_portfolio = []
    for _ in 1:n_pp_green
        green_pp = PowerPlant(
                    type = "Green",
                    c = 0.0,
                    freq = global_param.freq_per_powerplant,
                    capacity = global_param.freq_per_powerplant,
                    Aᵀ = 0.0,
                    em = 0.0
                   )
        push!(green_portfolio, green_pp)
    end

    # Initialize dirty power plants
    dirty_portfolio = []
    for _ in n_pp_green+1:n_pp
        dirty_pp = PowerPlant(
                    type = "Dirty",
                    c = 0.0,
                    freq = global_param.freq_per_powerplant,
                    capacity = global_param.freq_per_machine * init_param.Aᵀ_0,
                    Aᵀ = init_param.Aᵀ_0,
                    em = init_param.emᵀ_0
                   )
        update_c_pp!(dirty_pp, init_param.p_f)
        push!(dirty_portfolio, dirty_pp)
    end

    # Initialize ep struct
    energy_producer = EnergyProducer(
                        T = T,
                        green_portfolio = green_portfolio,
                        dirty_portfolio = dirty_portfolio,
                        Aᵀ_d = fill(init_param.Aᵀ_0, T),
                        emᵀ_d = fill(init_param.emᵀ_0, T)
                      )
    return energy_producer
end



"""
PRODUCTION
"""

"""
Production process of ep
"""
function produce_energy_ep!(
    ep::EnergyProducer,
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    t::Int,
    model::ABM
    )

    # Determine demand for energy
    compute_Dₑ_ep!(ep, all_cp, all_kp, t, model)

    # Check if production capacity needs to be expanded
    if ep.Dₑ[t] > ep.green_capacity + ep.dirty_capacity

    end

    choose_powerplants_ep!(ep, t)

    compute_Πₑ_ep!(ep, t)
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
    if ep.Dₑ[t] < ep.green_capacity
        ep.infra_marg = ep.green_portfolio
    else

        ep.infra_marg = ep.green_portfolio
        total_capacity = ep.green_capacity

        # Sort dirty portfolio as to take low cost power plants first
        sort!(ep.dirty_portfolio, by = pp -> pp.c)

        for dirty_pp in ep.dirty_portfolio
            push!(ep.infra_marg, dirty_pp)
            total_capacity += dirty_pp.capacity
            if total_capacity >= ep.Dₑ[t]
                break
            end
        end

    end
end


"""
INVESTMENT
"""


"""
Investment process of ep
"""
function plan_investments_ep!(
    ep::EnergyProducer, 
    t::Int
    )

    compute_Q̄_ep!(ep, t)

    compute_EIᵈ_ep!(ep, t)
    
end


"""

"""
function pick_expansion_investment_ep!(
    ep::EnergyProducer,
    bₑ::Float64,
    t::Int
    )

    if ep.IC_g[t] <= bₑ * ep.c_d[t]
        # Invest green
        ep.ECₑ[t] = ep.IC_g[t] * ep.EIᵈ[t]
        increase_n_powerplants_ep!(ep, type="Green")
    else
        # Invest dirty
        increase_n_powerplants_ep!(ep, type="Dirty")
    end

end


"""
Expands the amount of power plants based on decided expansion
"""
function increase_n_powerplants_ep!(
    ep::EnergyProducer;
    type::String
    )

    if type == "Green"
        # Check if current stock already has best type

            # Increase freq of best type

            # Create a new instance of best power plant

    else
        # Check if current stock already has best type

            # Increase freq of best type

            # Create a new instance of best power plant

    end
end


"""
INNOVATION
"""

"""
Innocation process
"""
function innovate_ep!(
    ep::EnergyProducer,
    global_param::GlobalParam,
    t::Int
    )

    # Compute R&D spending (Lamperti et al (2018), eq 18)
    ep.RDₑ[t] = t > 1 ? global_param.νₑ * ep.pₑ[t-1] * ep.Dₑ[t-1] : 0.0

    # Compute portions of R&D spending going to innovation in green and dirty tech
    #   (Lamperti et al (2018), eq 18.5)
    ep.IN_g[t] = global_param.ξₑ * ep.RDₑ[t]
    ep.IN_d[t] = (1 - global_param.ξₑ) * ep.RDₑ[t]

    # Define success probabilities of tech search (Lamperti et al (2018), eq 19)
    θ_g = 1 - exp(global_param.ζ_ge * ep.IN_g[t])
    θ_d = 1 - exp(global_param.ζ_de * ep.IN_d[t])

    # Candidate innovation for green tech
    if rand(Bernoulli(θ_g))
        # Draw from Beta distribution
        κ_g = rand(Beta(global_param.α1, global_param.β1))

        # Scale over supports
        κ_g = global_param.κ_lower + κ_g * (global_param.κ_upper - global_param.κ_lower)

        # Compute possible cost, replace all future if better
        #   (Lamperti et al (2018), eq 20)
        poss_IC_g = ep.IC_g[t] * κ_g
        ep.IC_g[t:end] = poss_IC_g < ep.IC_g[t] ? poss_IC_g : ep.IC_g[t]
    end

    # Candidate innovation for dirty tech
    if rand(Bernoulli(θ_d))
        # Draw from Beta distribution
        κ_d_A, κ_d_em = rand(Beta(global_param.α1, global_param.β1), 2)

        # Scale over supports
        κ_d_A = global_param.κ_lower + κ_d_A * (global_param.κ_upper - global_param.κ_lower)
        κ_d_em = global_param.κ_lower + κ_d_em * (global_param.κ_upper - global_param.κ_lower)

        # Compute possible thermal efficiency, replace all future if better
        #   (Lamperti et al (2018), eq 21)
        poss_A_d = ep.Aᵀ_d[t] * (1 + κ_d_A)
        ep.Aᵀ_d[t:end] = poss_A_d > ep.Aᵀ_d[t] ? poss_A_d : ep.Aᵀ_d[t]

        # Compute possible emissions, replace all future if better
        #   (Lamperti et al (2018), eq 21)
        poss_em_d = ep.emᵀ_d[t] * (1 - κ_d_em)
        ep.emᵀ_d[t:end] = poss_em_d < ep.emᵀ_d[t] ? poss_em_d : ep.emᵀ_d[t]
    end
end



"""
COMPUTE AND UPDATE FUNCTIONS
"""


"""
Computes the profits resulting from last periods production.
    Lamperti et al (2018), eq 10.
"""
function compute_Πₑ_ep!(
    ep::EnergyProducer, 
    t::Int
    )

    ep.Πₑ[t] = ep.pₑ[t] * ep.Dₑ[t] - ep.PCₑ[t] - ep.ICₑ[t] - ep.RDₑ[t]
end


"""
Computes cost of production PCₑ of infra marginal power plants
    Lamperti et al (2018), eq 12.
"""
function compute_PCₑ_ep!(
    ep::EnergyProducer,
    t::Int
    )

    dirty_cost = 0.0

    for pp in ep.infra_marg
        if pp ∈ ep.dirty_portfolio
            dirty_cost += pp.freq * pp.c * pp.Aᵗ
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

    c̄ = ep.dirty_portfolio[end].c
    ep.pₑ[t] = ep.Dₑ[t] <= ep.green_capacity ? ep.μₑ : c̄ + ep.μₑ
end


"""
Updates capacity figures for green and dirty technologies
    Lamperti et al (2018) eq 14.
"""
function update_capacities_ep!(
    ep::EnergyProducer
    )

    ep.green_capacity = sum(pp -> pp.freq, ep.green_portfolio)
    ep.dirty_capacity = sum(pp -> pp.freq, ep.dirty_portfolio)
end


"""
Computes the maximum production level Q̄
    Lamperti et al (2018) eq 15.
"""
function compute_Q̄_ep!(
    ep::EnergyProducer, 
    t::Int
    )

    ep.Q̄ₑ[t] = ep.green_capacity + sum(pp -> pp.freq * pp.Aᵀ, ep.dirty_portfolio)
end


"""
Computes the desired expansion investment
    Lamperti et al (2018) eq 16.
"""
function compute_EIᵈ_ep!(
    ep::EnergyProducer,
    t::Int
    )

    # TODO: I dont know yet what Kᵈ is, likely has to be changed
    ep.EIᵈ[t] = ep.Q̄ₑ[t] < ep.Dₑ[t] ? cp.Dₑ[t] - (cp.green_capacity + cp.dirty_capacity) : 0
end


"""
Computes total energy demand in period
"""
function compute_Dₑ_ep!(
    ep::EnergyProducer,
    all_cp::Vector{Int},
    all_kp::Vector{Int},
    t::Int,
    model::ABM
    )

    ep.Dₑ[t] = sum(cp_id -> model[cp_id].EU, all_cp) + sum(kp_id -> model[kp_id].EU, all_kp)
end


"""
Updates age of power plants
"""
function update_age_pp_ep!(
    ep::EnergyProducer
    )

    for pp in Iterators.flatten((ep.green_portfolio, ep.dirty_portfolio))
        pp.age += 1
    end
end