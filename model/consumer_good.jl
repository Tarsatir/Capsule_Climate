using Statistics


"""
Defines struct for consumer good producer
"""
mutable struct ConsumerGoodProducer <: AbstractAgent
    id :: Int                   # id
    cp_id :: Int                # cp id 
    p :: Array{Float64}         # hist prices
    c :: Array{Float64}         # hist cost
    RD :: Array{Float64}        # R&D spending
    D :: Array{Float64}         # hist demand
    Dᵉ :: Float64               # exp demand
    N :: Array{Float64}         # hist inventory
    Nᵈ :: Float64               # desired inventory
    Q :: Array{Float64}         # hist production
    Qᵉ :: Float64               # exp production
    I :: Array{Float64}         # hist investments
    Ξ :: Array{Machine}         # capital stock
    Emp:: Array{AbstractAgent}  # employees list
    L :: Float64                # labor units
    Lᵉ:: Float64                # exp labor force
    ΔLᵈ :: Float64              # desired change in labor force
    w :: Array{Float64}         # wage level
    brochures :: Array          # brochures from kp
    π :: Array{Float64}         # hist productivity
    f :: Float64                # hist market share
    μ :: Array{Float64}         # hist markup
    Π :: Array{Float64}         # hist profits
    cI :: Float64               # internal funds for investments
    ΔDeb :: Float64             # changes in debt level
    balance :: Balance          # balance sheet
end


function initialize_cp(id :: Int, cp_id :: Int, machine_struct, n_consrgood :: Int)
    cp = ConsumerGoodProducer(
        id,                     # global id
        cp_id,                  # cp id
        [],                     # p: hist prices
        [],                     # c: hist cost
        [],                     # RD: hist R&D spending
        [36e3],                  # D: hist demand
        37e3,                    # Dᵉ exp demand
        [36e3],                 # N: hist inventory
        30e3,                   # Nᵈ: desired inventory 
        [36e3],                # Q: hist production
        37e3,                  # Qᵉ: exp production
        [rand()],               # I: hist investments
        [machine_struct],       # Ξ: capital stock
        [],                     # Emp: employees
        0,                      # L: labor units in company
        90,                      # Lᵉ: exp labor force
        0,                      # ΔLᵈ: desired change in labor force
        [1.0],                  # w: wage level
        [],                     # brochures from kp
        [rand()],               # π: hist productivity
        1/n_consrgood,          # f: market share
        [0.05],                 # μ: hist markup
        [100000],               # Π: hist profits
        0,                      # cI: internal funds for investments
        0,                      # ΔDeb: changes in debt level
        Balance(               
                30e3,         # - N: inventory
                0.0,            # - K: capital
                0.0,            # - NW: liquid assets
                0.0,            # - Deb: debt
                0.0             # - EQ: equity
            )
    )
    return cp
end


"""
Plans production amounts for consumer good producer (short term)
    - updates ST expected demand
    - determines ST production goals
    - based on ST, set labor demand
"""
function plan_production_cp!(cp :: AbstractAgent, global_param)

    # update expected demand
    cp.Dᵉ = global_param.ωD * cp.D[end] + (1 - global_param.ωD) * cp.Dᵉ

    # determine desired short-term production
    Qˢ = cp.Dᵉ + cp.Nᵈ - cp.balance.N

    # print(Qˢ)

    # compute corresponding change in labor stock
    total_prod = sum(map(machine -> machine.A * machine.freq, cp.Ξ))
    # println("L:", cp.L)
    # println(Qˢ/total_prod)
    cp.ΔLᵈ = Qˢ/total_prod - cp.L

    # update markup μ
    μ = compute_μ_cp(cp, global_param.υ, global_param.μ1)
    push!(cp.μ, μ)

    # update cost of production c
    c = compute_c_cp(cp, Qˢ)
    push!(cp.c, c)

    # compute price
    p = (1 + μ) * c
    # println(p)
    push!(cp.p, p)
end


"""
Plans production amounts for consumer good producer (long term)
    - updates LT expected demand
    - updates LT labor supply 
    - determines LT production goals
    - based on LT, set investment amount
"""
function plan_investment_cp!(cp :: AbstractAgent, global_param, all_kp :: Array{AbstractAgent})

    # choose kp
    p_star, c_star, kp_choice, cop_star, Aᵈ = choose_producer_cp(cp, global_param.b, all_kp)

    # plan replacement investments
    RS = plan_replacement_cp(cp, global_param, p_star, c_star)

    # update LT demand
    cp.Qᵉ = global_param.ωQ * cp.Q[end] + (1 - global_param.ωQ) * cp.Qᵉ

    # update expected labor supply
    cp.Lᵉ = global_param.ωL * cp.L[end] + (1 - global_param.ωL) * cp.Lᵉ

    # compute desired capital stock expansion
    Kᵈ = (cp.Qᵉ / cp.Lᵉ - sum(map(x -> x.A * x.freq, cp.Ξ)))/Aᵈ

    EIᵈ = Kᵈ - sum(map(x -> x.freq, cp.Ξ))
    
    Iₜ = EIᵈ + sum(map(x -> x.freq, RS))

    println("I", Iₜ)

    # TODO does not check if funds are available
    if Iₜ > 0
        order_machines_cp!(kp_choice, cp, Iₜ)
    end

end


# """
# Computes capital stock K and average productivity π_t

# # Input
# - consumer good producer struct

# # Output
# - capital stock at time t
# - average productivity at time t
# """
# function compute_K_π!(cp)

#     # Dosi et al (2013) Eq. 20.5
#     K_t = sum(map(x -> x.freq, cp.Ξ))
#     push!(K_t, cp.Bal.K)

#     # Dosi et al (2013) Eq. 21.5
#     π_t = sum(map(x -> x.A * (x.freq / Kt), cp.Ξ))
#     push!(π_t, cp.π)

#     return Kt, π_t
# end

# Dosi et al (2013) Eq. 17, computes cost of production
cop(p_t, c_t, b) = p_t + b * c_t


# function compute_p_c_E!(cp, global_params, w_t, l_t)

#     c_t = w_t / cp.π[end]
#     p_t = (1 + global_param.μ1) * c_t
#     push!(c_t, cp.c)
#     push!(p_t, cp.p)

#     E_t = -global_params.ω1 * p_t -global_params.ω2 * l_t
#     push!(E_t, cp.E)

# end


# function set_production!(cp, E_bar, χ, C_t, r)

#     # update debt levels
#     Deb = cp.Bal.Deb
#     ΔDeb = cp.ΔDeb
#     push!(Deb + ΔDeb, cp.Bal.Deb)

#     # compute new market share
#     E_t = cp.E[end]
#     # Dosi et al (2013) Eq. 13
#     f_t = cp.f[end] * (1 + χ * (E_t - E_bar)/E_bar)
#     push!(f_t, cp.f)

#     # compute true demand
#     D_t = f_t * C_t
#     push!(D_t, cp.D)

#     # compute profits
#     p_t = cp.p[end]
#     c_t = cp.c[end]
#     Q_t = cp.Q[end]
#     Deb_t = cp.Bal.Deb[end]
#     NW_t_1 = cp.Bal.NW[end]
#     cI = cp.cI

#     # Dosi et al (2013), after Eq. 13
#     Π_t = p_t * D_t - c_t * Q_t - r * Deb_t
#     push!(Π_t, cp.Π)

#     # Dosi et al (2013), after Eq. 13
#     NW_t = NW_t_1 + Π_t - cI
#     push!(NW_t, cp.Bal.NW)

# end


# function plan_investment!(cp, global_param, model_struct)

#     # choose capital goods producer
#     p_star, c_star, chosen_producer, cop_star = choose_producer_cp(cp.brochures, global_param.b)

#     # compute capital stock K and average productivity
#     K_t, π_t = compute_productivity!(cp)
#     push!(π_t, cp.π)

#     # compute investment
#     EId = plan_expansion(cp, global_param, K_t)
#     RS = plan_replacement_cp(cp, global_param, p_star, c_star)
#     I_t = EId + sum(map(x -> x.freq, RS)) 
#     push!(I_t, cp.I)

#     # determine how liquid assets and debt changes

#     #  TODO: repay debts
#     cp.cI = minimum(It, cp.Bal.NW[end])
#     cp.ΔDeb = maximum(It - cp.Bal.NW[end], 0)
    
#     order_machines_cp!(kp, cp, I_t)
# end


# function plan_expansion(cp, global_param, K_t)
#     Nd = global_param.ι * cp.D[end]
#     Qd = cp.D[end] + Nd - cp.Bal.N[end]
#     push!(Qd, cp.Q)
    
#     K_d = Qd / global_param.cu

#     # TODO: check if this is the correct function
#     EId = K_d - K_t

#     return EId
# end


function plan_replacement_cp(cp :: AbstractAgent, global_param, p_star :: Float64, c_star :: Float64)
    RS = []

    for machine in cp.Ξ
        if (p_star/(machine.A - c_star) <= global_param.b) || machine.age >= global_param.η
            push!(RS, machine)
        end
    end

    return RS
end


function choose_producer_cp(cp :: AbstractAgent, b :: Int, all_kp :: Array{AbstractAgent})

    if (length(cp.brochures) == 0)
        # in case of no brochures, pick a random kp
        chosen_producer = sample(all_kp)
        brochure = chosen_producer.brochure
        p_star = brochure[2]
        c_star = brochure[3]
        A_star = brochure[4]
        cop_star = cop(p_star, c_star, b)
        
        return p_star, c_star, chosen_producer, cop_star, A_star
    end

    # take first producer as best, then compare to other producers
    chosen_producer = cp.brochures[1][1]
    p_star = cp.brochures[1][2]
    c_star = cp.brochures[1][3]
    A_star = cp.brochures[1][4]
    cop_star = cop(p_star, c_star, b)

    for brochure in cp.brochures

        p_h = brochure[2]
        c_h = brochure[3]
        A_h = brochure[4]
        potential_cop = cop(p_h, c_h, b)

        # if cheaper, choose cheaper option
        if potential_cop < cop_star
            chosen_producer = brochure[1]
            p_star = p_h
            c_star = c_h
            A_star = A_h
            cop_star = potential_cop
        end

    end

    return p_star, c_star, chosen_producer, cop_star, A_star
end


function compute_μ_cp(cp :: AbstractAgent, υ :: Float64, μ1 :: Float64)
    μ = μ1
    if (length(cp.f) > 2)
        μ = cp.μ[end] * (1 + υ * (cp.f[end] - cp.f[end-1])/cp.f[end-1])
    end
    return μ
end

function compute_c_cp(cp :: AbstractAgent, Qˢ :: Float64)
    Ā = sum(map(x -> x.A * x.freq, cp.Ξ))
    c = (cp.w[end] * Qˢ / Ā) / Qˢ
    return c
end


function order_machines_cp!(kp_choice :: AbstractAgent, cp :: AbstractAgent, Iₜ :: Float64)
    order = (cp, Iₜ)
    push!(kp_choice.orders, order)
end

function reset_brochures_cp!(cp :: AbstractAgent)
    cp.brochures = []
end


function receive_machines!(RS, cp, new_machines)

end