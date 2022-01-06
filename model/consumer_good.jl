using Statistics


mutable struct Balance
    # assets
    N :: Array{Float64} # inventories
    K :: Array{Float64} # capital
    NW :: Array{Float64} # liquid assets

    # liabilities
    Deb :: Array{Float64} # debt
    EQ :: Array{Float64} # equity
end

"""
Defines struct for consumer good producer
"""
mutable struct ConsumerGoodProducer <: AbstractAgent
    id :: Int
    p :: Array{Float64}
    c :: Array{Float64}
    RD :: Array{Float64}
    D :: Array{Float64}
    Q :: Array{Float64}
    I :: Array{Float64}
    Ξ :: Array
    brochures :: Array
    π :: Array{Float64}
    f :: Array{Float64}
    Π :: Array{Float64}
    cI :: Float64
    ΔDeb :: Float64
    Balance :: Balance
end


"""
Computes capital stock K and average productivity π_t

# Input
- consumer good producer struct

# Output
- capital stock at time t
- average productivity at time t
"""
function compute_K_π!(cp)

    # Dosi et al (2013) Eq. 20.5
    K_t = sum(map(x -> x.freq, cp.Ξ))
    push!(K_t, cp.Bal.K)

    # Dosi et al (2013) Eq. 21.5
    π_t = sum(map(x -> x.A * (x.freq / Kt), cp.Ξ))
    push!(π_t, cp.π)

    return Kt, π_t
end

# Dosi et al (2013) Eq. 17, computes cost of production
cop(p_t, c_t, b) = p_t + b * c_t


function compute_p_c_E!(cp, global_params, w_t, l_t)

    c_t = w_t / cp.π[end]
    p_t = (1 + global_param.μ1) * c_t
    push!(c_t, cp.c)
    push!(p_t, cp.p)

    E_t = -global_params.ω1 * p_t -global_params.ω2 * l_t
    push!(E_t, cp.E)

end


function set_production!(cp, E_bar, χ, C_t, r)

    # update debt levels
    Deb = cp.Bal.Deb
    ΔDeb = cp.ΔDeb
    push!(Deb + ΔDeb, cp.Bal.Deb)

    # compute new market share
    E_t = cp.E[end]
    # Dosi et al (2013) Eq. 13
    f_t = cp.f[end] * (1 + χ * (E_t - E_bar)/E_bar)
    push!(f_t, cp.f)

    # compute true demand
    D_t = f_t * C_t
    push!(D_t, cp.D)

    # compute profits
    p_t = cp.p[end]
    c_t = cp.c[end]
    Q_t = cp.Q[end]
    Deb_t = cp.Bal.Deb[end]
    NW_t_1 = cp.Bal.NW[end]
    cI = cp.cI

    # Dosi et al (2013), after Eq. 13
    Π_t = p_t * D_t - c_t * Q_t - r * Deb_t
    push!(Π_t, cp.Π)

    # Dosi et al (2013), after Eq. 13
    NW_t = NW_t_1 + Π_t - cI
    push!(NW_t, cp.Bal.NW)

end


function plan_investment!(cp, global_param, model_struct)

    # choose capital goods producer
    p_star, c_star, chosen_producer, cop_star = choose_producer(cp.brochures, global_param.b)

    # compute capital stock K and average productivity
    K_t, π_t = compute_productivity!(cp)
    push!(π_t, cp.π)

    # compute investment
    EId = plan_expansion(cp, global_param, K_t)
    RS = plan_replacement(cp, global_param, p_star, c_star)
    I_t = EId + sum(map(x -> x.freq, RS)) 
    push!(I_t, cp.I)

    # determine how liquid assets and debt changes

    #  TODO: repay debts
    cp.cI = minimum(It, cp.Bal.NW[end])
    cp.ΔDeb = maximum(It - cp.Bal.NW[end], 0)
    
    order_machines!(kp, cp, I_t)
end


function plan_expansion(cp, global_param, K_t)
    Nd = global_param.ι * cp.D[end]
    Qd = cp.D[end] + Nd - cp.Bal.N[end]
    push!(Qd, cp.Q)
    
    K_d = Qd / global_param.cu

    # TODO: check if this is the correct function
    EId = K_d - K_t

    return EId
end


function plan_replacement(cp, global_param, p_star, c_star)
    RS = Array{Float64}

    for machine in cp.Ξ
        if (p_star/(machine.A - c_star) <= global_param.b) || machine.age >= global_param.η
            push!(RS, machine)
        end
    end

    return RS
end


function choose_producer(brochures, b)

    # take first producer as best, then compare to other producers
    chosen_producer = brochures[1][1]
    p_star = brochures[1][2]
    c_star = brochures[1][3]
    cop_star = cop(p_star, c_star, b)

    for brochure in brochures # TODO let skip first

        p_h = brochure[2]
        c_h = brochure[3]
        potential_cop = cop(p_h, c_h, b)

        # if cheaper, choose cheaper option
        if potential_cop < cop_star
            chosen_producer = brochure[1]
            p_star = p_h
            c_star = c_h
            cop_star = potential_cop
        end

    end

    return p_star, c_star, chosen_producer, cop_star
end


function order_machines!(kp, cp, I_t)
    order = (cp, I_t)
    push!(order, kp.orders)
end

function reset_brochures!(cp)
    cp.brochures = []
end


function receive_machines!(RS, cp, new_machines)

end