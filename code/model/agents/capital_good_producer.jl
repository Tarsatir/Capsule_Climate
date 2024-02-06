@with_kw mutable struct CapitalGoodProducer <: AbstractAgent

    id::Int64                               # global id
    kp_i::Int64                             # kp index
    age::Int64 = 0                          # firm age
    
    # Technology and innovation
    A_LP::Float64 = 1.0                   # labor prod sold machines
    A_EE::Float64 = 1.0                   # energy efficiency of sold machines
    A_EF::Float64 = 1.0                   # environmental friendlines of sold machines
    B_LP::Float64 = 1.0                   # labor prod own production
    B_EE::Float64 = 1.0                   # energy efficiency of own production
    B_EF::Float64 = 1.0                   # environmental friendlines of own production
    RD::Float64 = 0.0                     # R&D expenditure
    IM::Float64 = 0.0                     # hist immitation expenditure
    IN::Float64 = 0.0                     # hist innovation expenditure

    # Price and cost data
    μ::Vector{Float64}                    # markup rates
    p::Vector{Float64} = fill(1+μ[end], 3) # hist price data
    c::Vector{Float64} = ones(Float64, 3) # hist cost data
    true_c::Float64 = 0.0                 # true unit cost
    
    # Employment
    employees::Vector{Int64} = Int64[]      # employees in company
    mean_skill::Float64 = 1.0             # mean skill level of employees
    L::Float64 = 0.0                      # labor units in company
    L_RD::Float64 = 0.0                   # labor units used for R&D
    Lᵈ::Float64 = L                       # desired labor units in company
    ΔLᵈ::Float64 = 0.0                    # desired change in labor force
    w̄::Vector{Float64}                    # wage level
    wᴼ::Float64 = w̄[end]                  # offered wage
    wᴼ_max::Float64 = 0.0                 # maximum offered wage

    Oᵉ::Float64 = 250.                    # total amount of machine units expected to be ordered
    O::Float64 = 0.                       # total amount of machines units ordered
    O_unmet::Float64 = 0.                 # total amount of machines ordered in first round that could not be made
    prod_cap::Int64 = 0                   # total production capacity
    prod_queue::Dict = Dict{Int64, Int64}() # production queue of machines
    Q::Vector{Float64} = zeros(Float64, 3)# production quantities
    D::Vector{Float64} = zeros(Float64, 3)# hist demand
    EU::Float64 = 0.                      # energy use in the last period

    HC::Vector{Int64} = []                  # hist clients
    Π::Vector{Float64} = zeros(Float64, 3)# hist profits
    Πᵀ::Vector{Float64} = zeros(Float64, 3)# hist profits after tax
    debt_installments::Vector{Float64}    # installments of debt repayments
    f::Vector{Float64}                    # market share
    orders::Dict = Dict{Int64, Int64}()   # orders
    balance::Balance                      # balance sheet
    curracc::FirmCurrentAccount = FirmCurrentAccount() # current account
    
    emissions::Float64 = 0.0              # carbon emissions in last period                    
end


"""
Initializes kp agent, default is the heterogeneous state, otherwise properties are given
    as optional arguments.
"""
function initialize_kp(
    id::Int64, 
    kp_i::Int64,
    n_captlgood::Int64,
    b::Int64;
    NW=1000,
    A_LP=1.0,
    A_EE=1.0,
    A_EF=1.0,
    B_LP=1.0,
    B_EE=1.0,
    B_EF=1.0,
    μ=0.2,
    w̄=1.0,
    f=1/n_captlgood,
    )

    kp = CapitalGoodProducer(
        id = id,                        
        kp_i = kp_i,                               
        A_LP = A_LP,
        A_EE = A_EE,
        A_EF = A_EF,                        
        B_LP = B_LP,
        B_EE = B_EE,
        B_EF = B_EF,
        debt_installments = zeros(Float64, b+1),                       
        μ = fill(μ, 3),
        w̄ = fill(w̄, 3), 
        wᴼ = w̄,
        f = fill(f, 3),
        balance = Balance(NW=NW, EQ=NW)
    )
    return kp
end


"""
Checks if innovation is performed, then calls appropriate functions
"""
function innovate_kp!(
    kp::CapitalGoodProducer, 
    government::Government,
    globalparam, 
    all_kp::Vector{Int64}, 
    kp_distance_matrix::Array{Float64},
    w̄::Float64,
    t::Int64,
    ep,
    model::ABM,
    )

    # Determine levels of R&D, and how to divide under IN and IM
    set_RD_kp!(kp, globalparam.ξ, globalparam.ν)
    tech_choices = [(kp.A_LP, kp.A_EE, kp.A_EF, kp.B_LP, kp.B_EE, kp.B_EF)]

    # Determine innovation of machines (Dosi et al (2010); eq. 4)
    θ_IN = 1 - exp(-globalparam.ζ * kp.IN)
    if rand(Bernoulli(θ_IN))
        A_LP_t_in = update_techparam_p(kp.A_LP, globalparam)
        A_EE_t_in = update_techparam_p(kp.A_EE, globalparam)
        A_EF_t_in = update_techparam_p(kp.A_EF, globalparam; is_EF=true)

        B_LP_t_in = update_techparam_p(kp.B_LP, globalparam)
        B_EE_t_in = update_techparam_p(kp.B_EE, globalparam)
        B_EF_t_in = update_techparam_p(kp.B_EF, globalparam; is_EF=true)

        push!(tech_choices, (A_LP_t_in, A_EE_t_in, A_EF_t_in, B_LP_t_in, B_EE_t_in, B_EF_t_in))
    end

    # Determine immitation of competitors
    θ_IM = 1 - exp(-globalparam.ζ * kp.IM)
    if rand(Bernoulli(θ_IM))
        imitated_tech = imitate_technology_kp(kp, all_kp, kp_distance_matrix, model)
        push!(tech_choices, imitated_tech)
    end

    choose_technology_kp!(kp, government, w̄, globalparam, tech_choices, t, ep)
end


"""
Lets kp choose technology
"""
function choose_technology_kp!(
    kp::CapitalGoodProducer,
    government::Government,
    w̄::Float64,
    globalparam::GlobalParam,
    tech_choices,
    t::Int64,
    ep
    )

    update_μ_p!(kp, globalparam.ϵ_μ, t)

    # Make choice between possible technologies
    if length(tech_choices) > 1
        #   Lamperti et al (2018), eq 1 and 2, augmented for tax rates
        c_h_cp = map(tech -> cop(w̄, tech[1], government.τᴱ, ep.p_ep[t], tech[2], government.τᶜ, tech[3]), tech_choices)
        c_h_kp = map(tech -> cop(kp.w̄[end], tech[4], government.τᴱ, ep.p_ep[t], tech[5], government.τᶜ,  tech[6]), tech_choices)
 
        p_h = map(c -> (1 + kp.μ[end]) * c, c_h_kp)
        r_h = c_h_cp .* globalparam.b .+ p_h
        idx = argmin(r_h)

        # Update tech parameters
        kp.A_LP = tech_choices[idx][1]
        kp.A_EE = tech_choices[idx][2]
        kp.A_EF = tech_choices[idx][3]
        kp.B_LP = tech_choices[idx][4]
        kp.B_EE = tech_choices[idx][5]
        kp.B_EF = tech_choices[idx][6]
    end
end


"""
Creates brochures and sends to potential clients.
"""
function send_brochures_kp!(
    kp::CapitalGoodProducer,
    all_cp::Vector{Int64}, 
    globalparam,
    model::ABM;
    n_hist_clients=50::Int64
    )

    # Update brochure
    update_brochure!(kp, model)

    # Send brochure to historical clients
    # NO LONGER NEEDED, CP REMEMBER KP
    # for cp_id in kp.HC
    #     push!(model[cp_id].brochures, kp.brochure)
    # end

    # Select new clients, send brochure
    NC_potential = setdiff(all_cp, kp.HC)

    if length(kp.HC) == 0
        n_choices = n_hist_clients
    else
        n_choices = round(Int64, globalparam.γ * length(kp.HC))
    end
    
    # Send brochures to new clients
    NC = sample(NC_potential, min(n_choices, length(NC_potential)); replace=false)
    for cp_id in NC
        add_kp_cp!(model[cp_id], kp.id)
    end 
end


"""
Initializes kp brochure and adds to model properties
"""
function init_brochure!(
    kp::CapitalGoodProducer,
    model::ABM
)

    brochure =  Dict(
                        :price => kp.p[end],
                        :A_LP => kp.A_LP,
                        :A_EE => kp.A_EE,
                        :A_EF => kp.A_EF
                    )

    model.kp_brochures[Symbol(kp.id)] = brochure
end


"""
Updates kp brochure in model properties
"""
function update_brochure!(
    kp::CapitalGoodProducer,
    model::ABM
)

    model.kp_brochures[Symbol(kp.id)][:price] = kp.p[end]
    model.kp_brochures[Symbol(kp.id)][:A_LP] = kp.A_LP
    model.kp_brochures[Symbol(kp.id)][:A_EE] = kp.A_EE
    model.kp_brochures[Symbol(kp.id)][:A_EF] = kp.A_EF
end


"""
Uses inverse distances as weights for choice competitor to immitate
"""
function imitate_technology_kp(
    kp::CapitalGoodProducer, 
    all_kp::Vector{Int64}, 
    kp_distance_matrix, 
    model::ABM
    )::Tuple{Float64, Float64, Float64, Float64, Float64, Float64}

    weights = map(x -> 1/x, kp_distance_matrix[kp.kp_i,:])
    idx = sample(all_kp, Weights(weights))
    
    return model[idx].A_LP, model[idx].A_EE, model[idx].A_EF, model[idx].B_LP, model[idx].B_EE, model[idx].B_EF
end


"""
Lets kp update unit costs
"""
function compute_c_kp!(
    kp::CapitalGoodProducer,
    government::Government,
    p_ep::Float64
    )

    c_t = cop(kp.w̄[end], kp.B_LP, government.τᴱ, p_ep, kp.B_EE, government.τᶜ, kp.B_EF)
    if kp.Q[end] > 0
        c_t += kp.RD / kp.Q[end]
    end
    shift_and_append!(kp.c, c_t)
end


"""
Lets kp set price
"""
function compute_p_kp!(
    kp::CapitalGoodProducer, 
    )

    shift_and_append!(kp.p, (1 + kp.μ[end]) * kp.c[end])
end


"""
Lets kp compute the production capacity
"""
function update_prod_cap_kp!(
    kp::CapitalGoodProducer,
    globalparam::GlobalParam
)

    kp.prod_cap = ceil(Int64, (kp.L * kp.B_LP) / globalparam.freq_per_machine)
    # println(kp.L, " ", kp.prod_cap)
end


"""
Determines the level of R&D, and how it is divided under innovation (IN) 
and immitation (IM). based on Dosi et al (2010)
"""
function set_RD_kp!(
    kp::CapitalGoodProducer, 
    ξ::Float64, 
    ν::Float64
    )

    # Determine total R&D spending at time t, 
    if kp.age > 1
        # Based on sales 
        #   (Dosi et al, 2010; eq. 3)
        prev_S = kp.p[end] * kp.Q[end]
        kp.RD = prev_S > 0 ? ν * prev_S : ν * max(kp.balance.NW, 0)
    else
        # Based on liquid assets, given sales are zero
        kp.RD = ν * max(kp.balance.NW, 0)
    end

    # Cap R&D to possible innovation based on labor available in company
    kp.RD = min(kp.RD, kp.L * kp.w̄[end])

    # Decide fractions innovation (IN) and immitation (IM), 
    #   (Dosi et al, 2010; eq. 3.5)
    kp.IN = ξ * kp.RD
    kp.IM = (1 - ξ) * kp.RD
end


function hascapacity(
    kp::CapitalGoodProducer
    )

    return sum(values(kp.orders)) <= kp.Q[end] + 1000
end


"""
    receive_order_kp!(kp::CapitalGoodProducer, cp_id::Int64)

Lets kp receive orders, adds client as historical clients if it is not yet.
"""
function receive_order_kp!(
    kp::CapitalGoodProducer,
    cp_id::Int64,
    order_size::Int64,
    freq_per_machine::Int64
    )

    kp.orders[cp_id] = order_size
    kp.O += order_size * freq_per_machine

    # If cp not in HC yet, add as a historical client
    if cp_id ∉ kp.HC
        push!(kp.HC, cp_id)
    end
end


"""
Based on received orders, sets labor demand to fulfill production.
"""
function plan_production_kp!(
    kp::CapitalGoodProducer,
    globalparam::GlobalParam,
    model::ABM
    )

    # Update average wage level
    update_w̄_p!(kp, model)
    
    # Determine total amount of capital units to produce and amount of labor to hire
    kp.O = sum(values(kp.orders)) * globalparam.freq_per_machine

    # Determine amount of labor to hire
    update_Lᵈ!(kp, globalparam.λ)

    # Update maximum offered wage
    update_wᴼ_max_kp!(kp)
end

# """
# Based on expected demand, sets labor demand to fulfill production.
# """
# function plan_production_kp!(
#     kp::CapitalGoodProducer,
#     globalparam::GlobalParam,
#     model::ABM
#     )

#     # Update average wage level
#     update_w̄_p!(kp, model)

#     # Update expected orders
#     update_Oᵉ_kp!(kp, globalparam.ω)
    
#     # Determine total amount of capital units to produce and amount of labor to hire
#     # kp.O = sum(values(kp.orders)) * globalparam.freq_per_machine

#     # Determine amount of labor to hire
#     update_Lᵈ!(kp, globalparam.λ)

#     # Update maximum offered wage
#     update_wᴼ_max_kp!(kp)
# end


"""
Update expected amount of orders
"""
function update_Oᵉ_kp!(
    kp::CapitalGoodProducer, 
    ω::Float64
)

    kp.Oᵉ = ω * kp.Oᵉ + (1 - ω) * (kp.O + kp.O_unmet)
    # kp.Oᵉ = ω * (kp.O * 1.1) + (1 - ω) * kp.Oᵉ
    # println("   ", kp.O / 25, " ", kp.O_unmet / 25, " ", kp.Oᵉ / 25)
    # kp.O = 0.
end


"""
Update desired level of labor
"""
function update_Lᵈ!(
    kp::CapitalGoodProducer, 
    λ::Float64
    )

    # kp.Lᵈ = λ * kp.L + (1 - λ) * (kp.Oᵉ / kp.B_LP + kp.RD / kp.w̄[end])
    kp.Lᵈ = λ * kp.L + (1 - λ) * (kp.O / kp.B_LP + kp.RD / kp.w̄[end])
    # kp.Lᵈ = kp.Oᵉ / kp.B_LP + kp.RD / kp.w̄[end]
    kp.ΔLᵈ = max(kp.Lᵈ - kp.L, -kp.L)
end


"""
Lets kp add goods to the production queue, based on available labor supply
"""
function produce_goods_kp!(
    kp::CapitalGoodProducer,
    ep,
    globalparam::GlobalParam,
    t::Int64
    )

    # Determine what the total demand is, regardless if it can be satisfied
    D = kp.O
    shift_and_append!(kp.D, D)

    # Determine how much labor is needed to produce a full machine
    req_L = globalparam.freq_per_machine / kp.B_LP

    kp.prod_queue = Dict{Int64,Int64}()

    # Check if production is constrained
    if kp.L * kp.B_LP >= kp.O

        # Enough labor available, perform full production
        for (cp_id, n_ordered) in kp.orders

            if n_ordered != 0
                kp.prod_queue[cp_id] = n_ordered
            end
        end

    else
        # Production constrained, determine amount of production possible
        # and randomly select which machines to produce
        n_poss_prod = floor(Int64, kp.L / req_L)

        for (cp_id, n_ordered) in kp.orders

            if n_poss_prod > n_ordered && n_ordered != 0
                # Full production possible
                kp.prod_queue[cp_id] = n_ordered
                n_poss_prod -= n_ordered
            elseif n_poss_prod > 0 && n_poss_prod < n_ordered && n_ordered != 0
                # Partial production possible
                kp.prod_queue[cp_id] = n_poss_prod
                n_poss_prod = 0
            end

            if n_poss_prod == 0
                break
            end
        end
    end

    # Append total production amount of capital units
    Q = sum(values(kp.prod_queue)) * globalparam.freq_per_machine
    shift_and_append!(kp.Q, Q)

    # Update energy use from production
    update_EU_TCE_kp!(kp, ep.p_ep[t])

    update_emissions_kp!(kp)

    # Empty order queue
    for cp_id in keys(kp.orders)
        kp.orders[cp_id] = 0
    end
end


"""
Updates the energy use (EU) and total cost of energy (TCE) of kp
"""
function update_EU_TCE_kp!(
    kp::CapitalGoodProducer, 
    p_ep::Float64
    )

    kp.EU = kp.Q[end] / kp.B_EE
    kp.curracc.TCE = p_ep * kp.EU
end


"""
Update amount of emissions that result from energy usage EU. 
"""
function update_emissions_kp!(
    kp::CapitalGoodProducer
    )

    kp.emissions = kp.EU * kp.B_EF
end


"""
Sends orders from production queue to cp.
"""
function send_ordered_machines_kp!(
    kp::CapitalGoodProducer,
    ep,
    globalparam::GlobalParam,
    t::Int64,
    model::ABM
    )

    for (cp_id, n_machines) in kp.prod_queue

        if n_machines > 0

            # Produce machines in production queue, send to cp
            machines = initialize_machine_stock(
                            globalparam.freq_per_machine, 
                            n_machines;
                            p = kp.p[end], 
                            A_LP = kp.A_LP,
                            A_EE = kp.A_EE,
                            A_EF = kp.A_EF
                        )
            Iₜ = n_machines * globalparam.freq_per_machine * kp.p[end]
            receive_machines_cp!(model[cp_id], ep, machines, Iₜ, t)
        end

        kp.prod_queue[cp_id] = 0
    end
    
    # Update sales
    kp.curracc.S = kp.Q[end] * kp.p[end]
end


"""
Resets order queue
"""
function reset_order_queue_kp!(
    kp::CapitalGoodProducer
    )

    kp.orders = []
end


"""
Lets kp select cp as historical clients
"""
function select_HC_kp!(
    kp::CapitalGoodProducer, 
    all_cp::Vector{Int64};
    n_hist_clients=10::Int64
    )

    kp.HC = sample(all_cp, n_hist_clients; replace=false)
end


"""
Lets kp update maximum offered wage
"""
function update_wᴼ_max_kp!(
    kp::CapitalGoodProducer
    )
    
    kp.wᴼ_max = kp.B_LP * kp.p[end] 
end


"""
Filters out historical clients if they went bankrupt
"""
function remove_bankrupt_HC_kp!(
    kp::CapitalGoodProducer,
    bankrupt_cp::Vector{Int64}
    )

    filter!(cp_id -> cp_id ∉ bankrupt_cp, kp.HC)
end


"""
Updates market share of all kp.
"""
function update_marketshare_kp!(
    all_kp::Vector{Int64},
    model::ABM
    )

    kp_market = sum(kp_id -> model[kp_id].D[end], all_kp)

    for kp_id in all_kp
        if kp_market == 0
            f = 1 / length(all_kp)
        else
            f = model[kp_id].D[end] / kp_market
        end
        shift_and_append!(model[kp_id].f, f)
    end
end


# """
# Updates the markup rate μ
# """
# function update_μ_kp!(
#     kp::CapitalGoodProducer
#     )

#     shift_and_append!(kp.μ, kp.μ[end])
# end


"""
Replaces bankrupt kp with new kp. Gives them a level of technology and expectations
    from another kp. 
"""
function replace_bankrupt_kp!(
    bankrupt_kp::Vector{Int64},
    bankrupt_kp_i::Vector{Int64},
    all_kp::Vector{Int64},
    globalparam::GlobalParam,
    indexfund::IndexFund,
    initparam::InitParam,
    macro_struct::MacroEconomy,
    t::Int64,
    model::ABM
    )

    # TODO: describe in model

    # Check if not all kp have gone bankrupt, in this case, 
    # kp with highest NW will not be removed
    if length(bankrupt_kp) == length(all_kp)
        kp_id_max_NW = all_kp[1]
        for i in 2:length(all_kp)
            if model[all_kp[i]].curracc.NW > model[kp_id_max_NW].curracc.NW
                kp_id_max_NW = all_kp[i]
            end
        end
        bankrupt_kp = bankrupt_kp[!kp_id_max_NW]
    end

    # Get all nonbankrupt kp
    nonbankrupt_kp = filter(kp_id -> kp_id ∉ bankrupt_kp, all_kp)

    # Get the technology frontier
    A_LP_max = maximum(kp_id -> model[kp_id].A_LP, nonbankrupt_kp)
    A_EE_max = maximum(kp_id -> model[kp_id].A_EE, nonbankrupt_kp)
    A_EF_min = minimum(kp_id -> model[kp_id].A_EF, nonbankrupt_kp)

    B_LP_max = maximum(kp_id -> model[kp_id].B_LP, nonbankrupt_kp)
    B_EE_max = maximum(kp_id -> model[kp_id].B_EE, nonbankrupt_kp)
    B_EF_min = minimum(kp_id -> model[kp_id].B_EF, nonbankrupt_kp)

    # Compute the average stock of liquid assets of non-bankrupt kp
    avg_NW = mean(kp_id -> model[kp_id].balance.NW, nonbankrupt_kp)
    NW_coefficients = rand(Uniform(globalparam.φ3, globalparam.φ4),
                           length(bankrupt_kp))

    # Compute share of investments that can be paid from the investment fund                       
    all_req_NW = sum(avg_NW .* NW_coefficients)
    frac_NW_if = decide_investments_if!(indexfund, all_req_NW, t)

    # Re-use id of bankrupted company
    for (i, (kp_id, kp_i)) in enumerate(zip(bankrupt_kp, bankrupt_kp_i))
        # Sample a producer of which to take over the technologies, proportional to the 
        # quality of the technology
        tech_coeff = (globalparam.φ5 + rand(Beta(globalparam.α2, globalparam.β2)) 
                                        * (globalparam.φ6 - globalparam.φ5))

        new_A_LP = max(A_LP_max * (1 + tech_coeff), initparam.A_LP_0)
        new_A_EE = max(A_EE_max * (1 + tech_coeff), initparam.A_LP_0)
        new_A_EF = min(A_EF_min * (1 - tech_coeff), initparam.A_LP_0)

        new_B_LP = max(B_LP_max * (1 + tech_coeff), initparam.B_LP_0)
        new_B_EE = max(B_EE_max * (1 + tech_coeff), initparam.B_LP_0)
        new_B_EF = min(B_EF_min * (1 - tech_coeff), initparam.B_LP_0)

        NW_stock = NW_coefficients[i] * avg_NW

        # Initialize new kp
        new_kp = initialize_kp(
            kp_id, 
            kp_i, 
            length(all_kp),
            globalparam.b;
            NW = NW_stock,
            A_LP = new_A_LP,
            A_EE = new_A_EE,
            A_EF = new_A_EF,
            B_LP = new_B_LP,
            B_EE = new_B_EE,
            B_EF = new_B_EF,
            μ = macro_struct.markup_kp[t],
            w̄ = macro_struct.w_avg[t],
            f = 0.0
        )

        # Reset brochure
        update_brochure!(new_kp, model)

        # Borrow the remaining funds
        borrow_funds_p!(new_kp, (1 - frac_NW_if) * NW_stock, globalparam.b)

        # Add agent to model
        add_agent!(new_kp, model)
    end
end


