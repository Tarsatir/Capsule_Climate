@with_kw mutable struct EnergyProducer <: AbstractAgent
    T::Int64=T                                    # Total number of iterations

    D_ep::Vector{Float64} = zeros(Float64, T)     # Demand for energy units over time
    Qmax_ep::Vector{Float64} = zeros(Float64, T)     # Maximum production of units over time


    # Prices, cost and investments
    markup_ep::Float64                            # Markup to determine price
    profit_ep::Vector{Float64} = zeros(Float64, T) # Profits over time
    NW_ep::Vector{Float64} = zeros(Float64, T)    # Stock of liquid assets over time
    p_ep::Vector{Float64} = zeros(Float64, T)     # Price of energy over time
    PC_ep::Vector{Float64} = zeros(Float64, T)    # Cost of generating Dₑ(t) units of energy over time
    FU::Vector{Float64} = zeros(Float64, T)       # Number of fuel units used for production
    IC_ep::Vector{Float64} = zeros(Float64, T)    # Expansion and replacement investments over time
    RD_ep::Vector{Float64} = zeros(Float64, T)    # R&D expenditure over time
    IN_g::Vector{Float64} = zeros(Float64, T)     # R&D spending allocated to innovation in green tech
    IN_d::Vector{Float64} = zeros(Float64, T)     # R&D spending allocated to innovation in dirty tech
    EId_ep::Vector{Float64} = zeros(Float64, T)   # Desired amount of units of expansionary investments
    RSd_ep::Vector{Float64} = zeros(Float64, T)   # Desired amount of units of replacement investments
    EC_ep::Vector{Float64} = zeros(Float64, T)    # Cost of spansionary investment
    carbontax::Vector{Float64} = zeros(Float64, T)# Paid carbon taxes

    # Technological parameters
    IC_g::Vector{Float64}                       # Fixed investment cost of the cheapest new green power plant
    A_therm_ep::Vector{Float64} = zeros(Float64, T)   # Thermal efficiency of new power plants
    emnew_ep::Vector{Float64} = zeros(Float64, T)  # Emissions of new power plants
    c_d::Vector{Float64} = zeros(Float64, T)    # Discounted production cost of the cheapest dirty plant 

    # Owned power plants
    green_portfolio::Vector{PowerPlant}         # Portfolio of all green powerplants
    green_capacity::Vector{Float64} = zeros(Float64, T) # capacity of set of green powerplants
    dirty_portfolio::Vector{PowerPlant}         # Portfolio of all dirty powerplants
    dirty_capacity::Vector{Float64} = zeros(Float64, T) # Capacity of set of dirty powerplants,
    green_frac_prod::Vector{Float64} = zeros(Float64, T) # Green fraction of total production
    infra_marg::Vector{PowerPlant} = PowerPlant[]   # Infra-marginal plants
    pp_tb_replaced::Vector{PowerPlant} = PowerPlant[]   # Power plants to be replaced
    pp_tb_retired::Vector{PowerPlant} = PowerPlant[]   # Power plants to be retired

    emissions::Vector{Float64} = zeros(Float64, T)
end


function initialize_energy_producer(
    T::Int64,
    initparam::InitParam,
    τᶜ::Float64,
    globalparam::GlobalParam
    )::EnergyProducer

    # Initialize power plants
    n_pp = ceil(Int64, initparam.n_powerplants_init / globalparam.freq_per_powerplant)
    n_pp_green = ceil(Int64, initparam.frac_green * n_pp)

    # Initialize green power plants
    green_portfolio = []
    for _ in 1:n_pp_green
        green_pp = PowerPlant(
                    type = "Green",
                    age = sample(0:globalparam.ηₑ),
                    c = 0.0,
                    freq = globalparam.freq_per_powerplant,
                    capacity = globalparam.freq_per_powerplant,
                    Aᵀ = 0.0,
                    em = 0.0
                   )
        push!(green_portfolio, green_pp)
    end

    # Initialize dirty power plants
    dirty_portfolio = []
    for _ in n_pp_green+1:n_pp
        dirty_pp = init_powerplant(
                    "Dirty",
                    sample(0:globalparam.ηₑ),
                    0.0,
                    initparam.Aᵀ_0,
                    initparam.emᵀ_0,
                    globalparam
                   )
        update_c_pp!(dirty_pp, globalparam.p_f, τᶜ)
        push!(dirty_portfolio, dirty_pp)
    end

    # Initialize ep struct
    energy_producer = EnergyProducer(
                        T = T,
                        green_portfolio = green_portfolio,
                        dirty_portfolio = dirty_portfolio,
                        markup_ep = initparam.markup_ep,
                        A_therm_ep = fill(initparam.Aᵀ_0, T),
                        emnew_ep = fill(initparam.emᵀ_0, T),
                        IC_g = fill(initparam.IC_g_0, T)
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
    government::Government,
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    globalparam::GlobalParam,
    indexfund::IndexFund,
    frac_green::Float64,
    t::Int64,
    t_warmup::Int64,
    model::ABM
    )

    # Update age for all pp and cost figure for dirty pp
    for pp ∈ Iterators.flatten((ep.green_portfolio, ep.dirty_portfolio))

        update_age_pp!(pp)

        if pp ∈ ep.dirty_portfolio
            update_c_pp!(pp, globalparam.p_f, government.τᶜ)
        end
    end

    # Update production capacities
    update_capacities_ep!(ep, t)
    compute_Q̄_ep!(ep, t)
    
    # Update cost frontier of ep
    compute_c_ep!(ep, government, globalparam.p_f, t)

    # Determine demand for energy
    compute_Dₑ_ep!(ep, all_cp, all_kp, t, model)

    # Check if production capacity needs to be expanded and old pp replaced
    plan_investments_ep!(ep, government, globalparam, frac_green, t, t_warmup)
    
    # Choose pp to use in production
    choose_powerplants_ep!(ep, t, globalparam)

    compute_FU_ICₑ_ep!(ep, globalparam.p_f, t)
    compute_emissions_ep!(ep, t)

    pay_carbontax_ep!(ep, government, t)

    # Compute profits
    compute_Πₑ_NWₑ_ep!(ep, t)
    pay_dividends_ep!(ep, indexfund, t)
end

# """
# Lets ep choose power plants to produce energy demand with
# """

# function choose_powerplants_ep!(
#     ep::EnergyProducer,
#     t::Int64
#     )
    
#     # Compute the total green capacity available based on a maximum of 80% of demand
#     green_capacity_total = floor(min(0.8 * ep.D_ep[t], ep.green_capacity[t]))
#     # print("D_ep: ", ep.D_ep[t], "\n")
#     # print("green_capacity[t]: ", ep.green_capacity[t], "\n")
#     # print("green_capacity_total: ", green_capacity_total, "\n")
    

    
#     if !isempty(ep.green_portfolio) && !isnan(green_capacity_total / ep.green_capacity[t])
#         max_green_plants = Int(floor(green_capacity_total / ep.green_capacity[t]))
#     else
#         max_green_plants = 0
#     end
    
    
#     # If the total demand at time t can be met using only green tech, 
#     # use only green tech and set fraction of production to 1.0
#     if ep.D_ep[t] < green_capacity_total
#         ep.infra_marg = copy(ep.green_portfolio)
#         ep.green_frac_prod[t] = 1.0
#     else
#         # If green tech is insufficient, use a combination of green and dirty tech
#         ep.infra_marg = copy(ep.green_portfolio[1:min(end, max_green_plants)])
#         total_capacity = Int(green_capacity_total)

#         # Sort the dirty portfolio by cost to take the lowest cost power plants first
#         sort!(ep.dirty_portfolio, by = pp -> pp.c)

#         # Loop over the dirty power plants in the sorted dirty portfolio
#         for dirty_pp in ep.dirty_portfolio
#             # Add the dirty power plant to the infrastructure margin list
#             push!(ep.infra_marg, dirty_pp)
#             # Add the capacity of the dirty power plant to the total capacity
#             # of the infrastructure margin
#             total_capacity += dirty_pp.capacity * dirty_pp.Aᵀ
            
#             # If the total capacity of the infrastructure margin exceeds the 
#             # demand at time t, exit the loop
#             if total_capacity >= ep.D_ep[t]
#                 break
#             end
#         end

#         # Compute the fraction of production that will come from green tech
#         # using the length of the green portfolio and the length of the
#         # infrastructure margin list
#         #ep.green_frac_prod[t] = length(ep.green_portfolio[1:max_green_plants]) / length(ep.infra_marg)
#         ep.green_frac_prod[t] = length(ep.green_portfolio) / length(ep.infra_marg)
#         print("ep.green_frac_prod[t]: ", ep.green_frac_prod[t], "\n")
#         ep.green_frac_prod[t] = length(ep.green_portfolio[1:min(end, max_green_plants)]) / length(ep.infra_marg)
#         print(ep.green_portfolio)
#         print("ep.green_frac_prod_new[t]: ", ep.green_frac_prod[t], "\n")


#     end
# end

# function choose_powerplants_ep!(
#     ep::EnergyProducer,
#     t::Int64
#     )

#     # Check if all production can be done using green tech, if not, compute cost
#     # of production using dirty tech
#     if ep.D_ep[t] < ep.green_capacity[t]
#         ep.infra_marg = copy(ep.green_portfolio)
#         ep.green_frac_prod[t] = 1.0
#     else

#         ep.infra_marg = copy(ep.green_portfolio)
#         total_capacity = ep.green_capacity[t]

#         # Sort dirty portfolio as to take low cost power plants first
#         sort!(ep.dirty_portfolio, by = pp -> pp.c)

#         for dirty_pp in ep.dirty_portfolio
#             push!(ep.infra_marg, dirty_pp)
#             total_capacity += dirty_pp.capacity * dirty_pp.Aᵀ
#             if total_capacity >= ep.D_ep[t]
#                 break
#             end
#         end

#         ep.green_frac_prod[t] = length(ep.green_portfolio) / length(ep.infra_marg)
#         print("ep.green_frac_prod[t]: ", ep.green_frac_prod[t], t)
#     end
# end


function choose_powerplants_ep!(
    ep::EnergyProducer,
    t::Int64,
    globalparam::GlobalParam
    )
    



 
    shortened_portfolio, total_capacity  = shorten_portfolio!(ep.green_portfolio, ep.D_ep[t], float(globalparam.green_limit ))
    #print("Total capacity green: ", total_capacity, "\n")
    #print("Total demand: ", ep.D_ep[t]*globalparam.green_limit, "\n")

    if globalparam.green_limit*ep.D_ep[t] > ep.green_capacity[t]
        ep.infra_marg = copy(ep.green_portfolio)
        
    else  
        ep.infra_marg = copy(shortened_portfolio)
    
    end
    # Sort dirty portfolio as to take low cost power plants first
    sort!(ep.dirty_portfolio, by = pp -> pp.c)

    for dirty_pp in ep.dirty_portfolio
        push!(ep.infra_marg, dirty_pp)
        total_capacity += dirty_pp.capacity * dirty_pp.Aᵀ
        if total_capacity >= ep.D_ep[t]
            break
        end
    end
    #print the length of shortened_portfolio and green_portfolio
    #print("length of shortened_portfolio: ", length(shortened_portfolio), "\n")
    #print("length of green_portfolio: ", length(ep.green_portfolio), "\n")
    #print("length of infra_marg: ", length(ep.infra_marg), "\n")
    #print("length of dirty_portfolio: ", length(ep.dirty_portfolio), "\n")
    ep.green_frac_prod[t] = length(shortened_portfolio) / length(ep.infra_marg)
    #print ep.infra_marg
    #print("infra_marg", ep.infra_marg)
    #print("ep.green_frac_prod[t]: ", ep.green_frac_prod[t])
    
end

function shorten_portfolio!(
    green_portfolio::Vector{PowerPlant}, 
    D_ep::Float64, 
    green_limit::Float64)

    if isempty(green_portfolio)
        return green_portfolio, 0.0
    end

    sorted_portfolio = sort(green_portfolio, by = pp -> pp.capacity)
    total_capacity = sum(pp.capacity for pp in sorted_portfolio)
    green_capacity_limit = green_limit * D_ep
 
    
    shortened_portfolio = PowerPlant[]
    current_capacity = 0.0
    
    for pp in sorted_portfolio
        current_capacity += pp.capacity
        if current_capacity > green_capacity_limit
            break
        end
        push!(shortened_portfolio, pp)
    end
    
    return shortened_portfolio, current_capacity
end





function pay_dividends_ep!(
    ep::EnergyProducer, 
    indexfund::IndexFund, 
    t::Int64
    )

    # TODO: describe

    b = 3

    # ep should have at least enough NW to pay cost of production for b months plus an
    # investment in a green plant
    req_NW = b * ep.PC_ep[t] + ep.IC_g[t]

    # Pay expenses to if, as 'indeterminate' producer receives these expenses,
    dividends = ep.PC_ep[t] + ep.IC_ep[t] + ep.RD_ep[t]
    
    # Gather excess liquid assets
    dividends += max(ep.NW_ep[t] - req_NW, 0)

    ep.NW_ep[t] = min(ep.NW_ep[t], req_NW)

    receive_dividends_if!(indexfund, dividends)
end


"""
INVESTMENT
"""


"""
Investment process of ep
"""
function plan_investments_ep!(
    ep::EnergyProducer,
    government::Government,
    globalparam::GlobalParam,
    frac_green::Float64, 
    t::Int64,
    t_warmup::Int64
    )
    
    compute_RSᵈ_ep!(ep, globalparam.ηₑ, t, globalparam)

    compute_EIᵈ_ep!(ep, t)

    expand_and_replace_pp_ep!(ep, globalparam, frac_green, t, t_warmup)
end


"""

"""
function expand_and_replace_pp_ep!(
    ep::EnergyProducer,
    globalparam::GlobalParam,
    frac_green::Float64,
    t::Int64,
    t_warmup::Int64
    )
    
    n_add_pp = ceil(Int64, (ep.EId_ep[t] + ep.RSd_ep[t]) / globalparam.freq_per_powerplant)
    green_fraction = ep.green_capacity[t]/(ep.green_capacity[t]+ep.dirty_capacity[t])
    # Determine what share of additional powerplants should be green and dirty
    if t >= t_warmup
        # If after warmup period, ep can plan investment based on cheapest tech
        n_add_green_pp = 0
        n_add_dirty_pp = 0
        
        if  ep.IC_g[t] <= globalparam.bₑ * ep.c_d[t] && green_fraction < globalparam.green_limit # && (ep.green_capacity[t]/(ep.green_capacity[t]+ep.dirty_capacity[t])) <= globalparam.green_limit
            #(floor(exp(green_fraction*10)) + ep.IC_g[t]) <= globalparam.bₑ * ep.c_d[t] &&
            #variable = (floor(exp(green_fraction*10)) + ep.IC_g[t])
            #print('\n' ," Check Condition", variable, "<=", globalparam.bₑ * ep.c_d[t], "\n")
            #print('\n' ," Check Condition", green_fraction, "<=", globalparam.green_limit, "\n")
            n_add_green_pp = n_add_pp
        else
            n_add_dirty_pp = n_add_pp
        end
    else
        # If in warmup period, ep will replace such that share of green tech remains
        # at init level

        # n_add_green_pp = max(ceil(frac_green * (ep.green_capacity[t] + ep.dirty_capacity[t]
        #                       + n_add_pp) - ep.green_capacity[t]), 0)
        # n_add_dirty_pp = max(n_add_pp - n_add_green_pp, 0)
        n_add_green_pp = round(Int64, frac_green * n_add_pp)
        n_add_dirty_pp = round(Int64, (1-frac_green) * n_add_pp)
    end

    # println(ep.EId_ep[t], " ", ep.RSd_ep[t])
    # println(n_add_pp, " ", n_add_green_pp, " ", n_add_dirty_pp)
    # println()

    # TODO: if expanding green, invested sum should go somewhere!

    # Add green power plants
    if n_add_green_pp > 1
        
        # Invest green
        #ep.EC_ep[t] = ep.IC_g[t] * ep.EId_ep[t] #FIXME
        ep.EC_ep[t] = (floor(exp(ep.green_frac_prod[t]*10))  + ep.IC_g[t]) * ep.EId_ep[t]

        
        # Build new green pp
        for _ in 1:n_add_green_pp
            green_pp = init_powerplant(
                "Green",
                0,
                0.,
                0.,
                0.,
                globalparam
            )
            # green_pp = PowerPlant(
            #     type = "Green",
            #     age = 0,
            #     c = 0.0,
            #     freq = globalparam.freq_per_powerplant,
            #     capacity = globalparam.freq_per_powerplant,
            #     Aᵀ = 0.0,
            #     em = 0.0
            # )
            push!(ep.green_portfolio, green_pp)
        end
    end

    # Add dirty power plants
    if n_add_dirty_pp > 1

        # Invest dirty
        for _ in 1:n_add_dirty_pp
            dirty_pp = init_powerplant(
                            "Dirty",
                            0,
                            ep.c_d[t],
                            ep.A_therm_ep[t],
                            ep.emnew_ep[t],
                            globalparam
                       )
            push!(ep.dirty_portfolio, dirty_pp)
        end
    end

    # Discard old power plants
    filter!(pp -> pp ∉ ep.pp_tb_replaced, ep.green_portfolio)
    filter!(pp -> pp ∉ ep.pp_tb_retired, ep.green_portfolio)

    filter!(pp -> pp ∉ ep.pp_tb_replaced, ep.dirty_portfolio)
    filter!(pp -> pp ∉ ep.pp_tb_retired, ep.dirty_portfolio)
end


"""
INNOVATION
"""

"""
Innocation process
"""
# function innovate_ep!(
#     ep::EnergyProducer,
#     globalparam::GlobalParam,
#     t::Int64
#     )

#     # Compute R&D spending (Lamperti et al (2018), eq 18)
#     ep.RD_ep[t] = t > 1 ? globalparam.νₑ * ep.p_ep[t-1] * ep.D_ep[t-1] : 0.0

#     # Compute portions of R&D spending going to innovation in green and dirty tech
#     #   (Lamperti et al (2018), eq 18.5)
#     ep.IN_g[t] = globalparam.ξₑ * ep.RD_ep[t]
#     ep.IN_d[t] = (1 - globalparam.ξₑ) * ep.RD_ep[t]

#     # Define success probabilities of tech search (Lamperti et al (2018), eq 19)
#     θ_g = 1 - exp(-globalparam.ζ_ge * ep.IN_g[t])
#     θ_d = 1 - exp(-globalparam.ζ_de * ep.IN_d[t])

#     # Candidate innovation for green tech
#     if rand(Bernoulli(θ_g))
#         # Draw from Beta distribution
#         κ_g = rand(Beta(globalparam.α1, globalparam.β1))

#         # Scale over supports
#         κ_g = globalparam.κ_lower + κ_g * (globalparam.κ_upper - globalparam.κ_lower)

#         # Compute possible cost, replace all future if better
#         #   (Lamperti et al (2018), eq 20)
#         poss_IC_g = ep.IC_g[t] * (1 - κ_g)
#         ep.IC_g[t:end] .= poss_IC_g < ep.IC_g[t] ? poss_IC_g : ep.IC_g[t]
#     end

#     # Candidate innovation for dirty tech
#     if rand(Bernoulli(θ_d))
#         # Draw from Beta distribution
#         κ_d_A, κ_d_em = rand(Beta(globalparam.α1, globalparam.β1), 2)

#         # Scale over supports
#         κ_d_A = globalparam.κ_lower + κ_d_A * (globalparam.κ_upper - globalparam.κ_lower)
#         κ_d_em = globalparam.κ_lower + κ_d_em * (globalparam.κ_upper - globalparam.κ_lower)

#         # Compute possible thermal efficiency, replace all future if better
#         #   (Lamperti et al (2018), eq 21)
#         poss_A_d = ep.A_therm_ep[t] * (1 + κ_d_A)
#         ep.A_therm_ep[t:end] .= poss_A_d > ep.A_therm_ep[t] ? poss_A_d : ep.A_therm_ep[t]

#         # Compute possible emissions, replace all future if better
#         #   (Lamperti et al (2018), eq 21)
#         poss_em_d = ep.emnew_ep[t] * (1 - κ_d_em)
#         ep.emnew_ep[t:end] .= poss_em_d < ep.emnew_ep[t] ? poss_em_d : ep.emnew_ep[t]
#     end
# end

function innovate_ep!(
    ep::EnergyProducer,
    globalparam::GlobalParam,
    t::Int64
    )

    # Compute R&D spending (Lamperti et al (2018), eq 18)
    ep.RD_ep[t] = t > 1 ? globalparam.νₑ * ep.p_ep[t-1] * ep.D_ep[t-1] : 0.0

    # Compute portions of R&D spending going to innovation in green and dirty tech
    #   (Lamperti et al (2018), eq 18.5)
    ep.IN_g[t] = globalparam.ξₑ * ep.RD_ep[t]
    ep.IN_d[t] = (1 - globalparam.ξₑ) * ep.RD_ep[t]

    # Define success probabilities of tech search (Lamperti et al (2018), eq 19)
    θ_g = 1 - exp(-globalparam.ζ_ge * ep.IN_g[t])
    θ_d = 1 - exp(-globalparam.ζ_de * ep.IN_d[t])

    # Check if θ_g is valid, otherwise opt for trivial innovation
    if isnan(θ_g) || θ_g < 0 || θ_g > 1
        println("Warning: Invalid θ_g = $θ_g for green tech innovation, opting for trivial innovation.")
        θ_g = 0  # Trivial innovation: skip the innovation process
    end

    # Candidate innovation for green tech
    if rand(Bernoulli(θ_g))
        # Draw from Beta distribution
        κ_g = rand(Beta(globalparam.α1, globalparam.β1))

        # Scale over supports
        κ_g = globalparam.κ_lower + κ_g * (globalparam.κ_upper - globalparam.κ_lower)

        # Compute possible cost, replace all future if better
        #   (Lamperti et al (2018), eq 20)
        poss_IC_g = ep.IC_g[t] * (1 - κ_g)
        ep.IC_g[t:end] .= poss_IC_g < ep.IC_g[t] ? poss_IC_g : ep.IC_g[t]
    end

    # Check if θ_d is valid, otherwise opt for trivial innovation
    if isnan(θ_d) || θ_d < 0 || θ_d > 1
        println("Warning: Invalid θ_d = $θ_d for dirty tech innovation, opting for trivial innovation.")
        θ_d = 0  # Trivial innovation: skip the innovation process
    end

    # Candidate innovation for dirty tech
    if rand(Bernoulli(θ_d))
        # Draw from Beta distribution
        κ_d_A, κ_d_em = rand(Beta(globalparam.α1, globalparam.β1), 2)

        # Scale over supports
        κ_d_A = globalparam.κ_lower + κ_d_A * (globalparam.κ_upper - globalparam.κ_lower)
        κ_d_em = globalparam.κ_lower + κ_d_em * (globalparam.κ_upper - globalparam.κ_lower)

        # Compute possible thermal efficiency, replace all future if better
        #   (Lamperti et al (2018), eq 21)
        poss_A_d = ep.A_therm_ep[t] * (1 + κ_d_A)
        ep.A_therm_ep[t:end] .= poss_A_d > ep.A_therm_ep[t] ? poss_A_d : ep.A_therm_ep[t]

        # Compute possible emissions, replace all future if better
        #   (Lamperti et al (2018), eq 21)
        poss_em_d = ep.emnew_ep[t] * (1 - κ_d_em)
        ep.emnew_ep[t:end] .= poss_em_d < ep.emnew_ep[t] ? poss_em_d : ep.emnew_ep[t]
    end
end




"""
COMPUTE AND UPDATE FUNCTIONS
"""


"""
Computes the profits and new liquid assets resulting from last periods production.
    Lamperti et al (2018), eq 10.
"""
function compute_Πₑ_NWₑ_ep!(
    ep::EnergyProducer, 
    t::Int64
    )

    # Compute profits
    ep.profit_ep[t] = ep.p_ep[t] * ep.D_ep[t] - ep.PC_ep[t] - ep.IC_ep[t] - ep.RD_ep[t] - ep.carbontax[t]
    ep.NW_ep[t] = t > 1 ? ep.NW_ep[t-1] + ep.profit_ep[t] : ep.profit_ep[t]
end


"""
Computes cost of production PCₑ of infra marginal power plants
    Lamperti et al (2018), eq 12.
"""
function compute_PCₑ_ep!(
    ep::EnergyProducer,
    t::Int64
    )

    dirty_cost = 0.0

    for pp in ep.infra_marg
        if pp ∈ ep.dirty_portfolio
            dirty_cost += pp.freq * pp.c * pp.Aᵗ
        end
    end

    ep.PC_ep[t] = dirty_cost
end


"""
Computes the price for each sold energy unit
"""
function compute_pₑ_ep!(
    ep::EnergyProducer, 
    t::Int64
    )

    # c̄ = length(ep.dirty_portfolio) > 0 ? ep.dirty_portfolio[end].c : 0.0
    # ep.p_ep[t] = t > 1 && ep.D_ep[t-1] <= ep.green_capacity[t] ? ep.markup_ep : c̄ + ep.markup_ep

    if length(ep.dirty_portfolio) > 0
        c̄ = maximum(pp -> pp.c, ep.dirty_portfolio)
    else
        c̄ = 0.
    end

    ep.p_ep[t] = t > 1 ? ep.markup_ep + (1 - ep.green_frac_prod[t-1]) * c̄ : ep.markup_ep
end


"""
Updates capacity figures for green and dirty technologies
    Lamperti et al (2018) eq 14.
"""
function update_capacities_ep!(
    ep::EnergyProducer,
    t::Int64
    )

    ep.green_capacity[t] = length(ep.green_portfolio) > 0 ? sum(pp -> pp.freq, ep.green_portfolio) : 0.0
    ep.dirty_capacity[t] = length(ep.dirty_portfolio) > 0 ? sum(pp -> pp.freq, ep.dirty_portfolio) : 0.0
end


"""
Computes the maximum production level Q̄
    Lamperti et al (2018) eq 15.
"""
function compute_Q̄_ep!(
    ep::EnergyProducer, 
    t::Int64
    )

    ep.Qmax_ep[t] = ep.green_capacity[t] 

    if length(ep.dirty_portfolio) > 0
        ep.Qmax_ep[t] += length(ep.dirty_portfolio) > 0 ? sum(pp -> pp.freq * pp.Aᵀ, ep.dirty_portfolio) : 0.0
    end
end


"""
Computes the desired expansion investment
    Lamperti et al (2018) eq 16.
"""
function compute_EIᵈ_ep!(
    ep::EnergyProducer,
    t::Int64
    )

    ep.EId_ep[t] = ep.Qmax_ep[t] < ep.D_ep[t] ? ep.D_ep[t] - (ep.green_capacity[t] + ep.dirty_capacity[t]) : 0.0
end


"""
Computes the desired replacement investment

"""
function compute_RSᵈ_ep!(
    ep::EnergyProducer,
    ηₑ::Int64,
    t::Int64,
    globalparam::GlobalParam
    )
    
    # Select powerplants to be replaced
    ep.pp_tb_replaced = []
    ep.pp_tb_retired = []

   
    dirty_excess_capacity =  ep.green_capacity[t] - globalparam.green_limit*ep.D_ep[t] #ep.D_ep[t]*(1-globalparam.green_limit)   # /ep.Qmax_ep[t]
    green_excess_capacity = ep.dirty_capacity[t] - (1-globalparam.green_limit)*ep.D_ep[t]#ep.D_ep[t]*globalparam.green_limit
    
    # Loop over all power plants (both green and dirty)
    for pp in Iterators.flatten((ep.green_portfolio, ep.dirty_portfolio))
        # Check if the age of the power plant is greater than or equal to the age limit
        if pp.age >= ηₑ
            # Choose the correct excess capacity depending on the type of power plant
            excess_capacity = pp.type == "Green" ? green_excess_capacity : dirty_excess_capacity

            # If the excess capacity is greater than the capacity of the power plant
            if excess_capacity > pp.capacity
                # Subtract the capacity of the power plant from the excess capacity
                excess_capacity -= pp.capacity
                # Add the power plant to the list of power plants to be retired
                push!(ep.pp_tb_retired, pp)
            else
                # If the excess capacity is not greater than the capacity of the power plant,
                # add the power plant to the list of power plants to be replaced
                push!(ep.pp_tb_replaced, pp)
            end
        end
    end

    # If there are any power plants to be replaced, calculate the total capacity
    # of those power plants. If not, set the desired replacement to 0.
    ep.RSd_ep[t] = length(ep.pp_tb_replaced) > 0 ? sum(pp -> pp.capacity, ep.pp_tb_replaced) : 0.0
end
#     for pp in Iterators.flatten((ep.green_portfolio, ep.dirty_portfolio))
#         if pp.age >= ηₑ
#             if excess_capacity > pp.capacity
#                 excess_capacity -= pp.capacity
#                 push!(ep.pp_tb_retired, pp)
#             else
#                 push!(ep.pp_tb_replaced, pp)
#             end
#         end
#     end
#     ep.RSd_ep[t] = length(ep.pp_tb_replaced) > 0 ? sum(pp -> pp.capacity, ep.pp_tb_replaced) : 0.0
# end


"""
Computes total energy demand in period
"""
function compute_Dₑ_ep!(
    ep::EnergyProducer,
    all_cp::Vector{Int64},
    all_kp::Vector{Int64},
    t::Int64,
    model::ABM
    )

    ep.D_ep[t] = sum(cp_id -> model[cp_id].EU, all_cp) + sum(kp_id -> model[kp_id].EU, all_kp)
end


"""
Computes cost of dirty powerplant at technological frontier
"""
function compute_c_ep!(
    ep::EnergyProducer,
    government::Government, 
    p_f::Float64, 
    t::Int64
    )

    ep.c_d[t] = p_f / ep.A_therm_ep[t] + government.τᶜ * ep.emnew_ep[t]
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

"""
Computes emissions following from energy production
"""
function compute_emissions_ep!(
    ep::EnergyProducer, 
    t::Int64
    )
   
    req_capacity = ep.D_ep[t]*(1-ep.green_frac_prod[t])
    #print("req_capacity: ", req_capacity, "\n")
    total_emissions = 0.0

    # Only use fraction of machines required, in order as sorted for costs
    #print("infra.marg", length(ep.infra_marg),"\n")
    for pp in ep.infra_marg
        if pp ∈ ep.dirty_portfolio
            total_emissions += max(min(req_capacity / pp.capacity, 1.0), 0) * pp.capacity * pp.em
            req_capacity -= pp.capacity
        end
    end

    ep.emissions[t] = total_emissions
    #print("Emissions: ", ep.emissions[t], "\n")
    
end


function pay_carbontax_ep!(
    ep::EnergyProducer,
    government::Government,
    t::Int64
    )

    ep.carbontax[t] = government.τᶜ * ep.emissions[t]
    # TODO move this to a gov function
    government.rev_carbontax[t] += ep.carbontax[t]
end


function compute_FU_ICₑ_ep!(
    ep::EnergyProducer,
    p_f::Float64, 
    t::Int64
)

    ep.FU[t] = length(ep.infra_marg) > 0 ? sum(pp-> pp ∈ ep.dirty_portfolio ? pp.capacity / pp.Aᵀ : 0.0, ep.infra_marg) : 0.0
    ep.IC_ep[t] = p_f * ep.FU[t]
end