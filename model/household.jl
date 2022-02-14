mutable struct Household <: AbstractAgent
    id :: Int                   # global id
    employed :: Bool            # is employed
    employer :: Int             # id of employer
    I :: Vector{Float64}        # hist income
    Iᵀ :: Vector{Float64}       # hist taxed income
    Iᵉ :: Float64               # expected income
    L :: Float64                # labor units in household
    S :: Array{Float64}         # total savings
    Sᵈ :: Float64               # desired savings
    s :: Float64                # savings rate
    B :: Array{Float64}         # budget
    C :: Float64                # cash
    N_B_min :: Float64          # substistence level of basic goods
    w :: Vector{Float64}        # wage
    wˢ :: Float64               # satisfying wage
    wʳ :: Float64               # requested wage
    wᵉ :: Float64               # expected wage
    ωI :: Float64               # memory param income expectation
    pref_bp_id                  # prefered basic good provider
    pref_lp_id                  # prefered luxury good provider
    T_unemp :: Int              # time periods unemployed
end

function initialize_hh(
    id::Int,
    τᴵ::Float64
    )::Household
    hh = Household(
        id,                     # global id
        false,                  # bool: employed
        0,                      # id of employer
        [],                     # I: hist income
        [],                     # Iᵀ: hist taxed income
        100*(1-τᴵ),             # Iᵉ: exp income
        100,                    # L: labor units
        [10],                   # S: total savings
        0,                      # Sᵈ: desired savings
        0,                      # s: savings rate
        [0.0],                  # B: budget
        100,                    # C: cash
        10,                     # N_B_min: substistence level of basic goods
        [1.0],                  # w: wage
        1.0,                    # wˢ: satisfying wage
        1.0,                    # wʳ: requested wage
        1.0,                    # wᵉ: expected wage
        0.5,                    # ωI: memory param income exp
        nothing,                # bg: prefered basic good provider
        nothing,                # lg: prefered luxury good provider
        0,                      # T_unemp: time periods unemployed
    )
    return hh
end


function pick_cp_hh!(
    hh::Household, 
    supplying_bp::Vector{Int}, 
    supplying_lp::Vector{Int}
    )

    # println(length(supplying_bp), " ", length(supplying_lp))
    
    # TODO: do this in a more sophisticated way
    hh.pref_bp_id = sample(supplying_bp)
    hh.pref_lp_id = sample(supplying_lp)

end


"""
Computes the expected income based on perceived probabilities.
"""
function compute_exp_income_hh!(
    hh::Household, 
    P_HU::Float64, 
    P_UU::Float64, 
    UB::Float64,
    model::ABM
    )

    # determine income expectation for employed workers
    if hh.employed

        # update expected wage
        # TODO: find a way to reset this when unemployed
        if length(hh.w) > 1
            hh.wᵉ = hh.ωI * hh.wᵉ + (1-hh.ωI) * (2*hh.w[end] - hh.w[end-1])
        end

        P_UE = model[hh.employer].P_FE * (1 - P_HU)

        hh.Iᵉ = P_UE * UB + (1 - P_UE) * hh.wᵉ

    # determine income expectation for unemployed workers
    else

        hh.Iᵉ = P_UU * UB + (1 - P_UU) * hh.wʳ

    end
end


"""
Determines savings rate s
"""
function set_savingsrate_hh!(
    hh::Household, 
    avg_T_unemp::Float64, 
    UB::Float64
    )

    # determine average budget over last year or over available information
    if length(hh.B[end]) >= 5
        B̄ = mean(hh.B[end-4, end])
    else
        B̄ = mean(hh.B)
    end
    
    if hh.employed
        # determine desired level of savings
        hh.Sᵈ = avg_T_unemp * (B̄ - UB)

        # determine savings rate
        s = max((hh.Sᵈ - hh.S[end]) / (hh.I[end] + 3*hh.Iᵉ), -hh.S[end]/hh.I[end])
    else

        s = (B̄ - UB) / UB

    end

    # TODO: find a solution for this
    if isnan(s)
        hh.s = 0
    else
        hh.s = s
    end

    # println(hh.s, " ", hh.I[end], " ", 3*hh.Iᵉ)

end


"""
Determines the optimal consumption package
    - Choose which goods to buy
    - Choose optimal consumption package for both goods
"""
function set_cons_package_hh!(
    hh::AbstractAgent, 
    τˢ::Float64,
    model::ABM
    )::Tuple{Float64, Float64}

    pref_bp = model[hh.pref_bp_id]
    pref_lp = model[hh.pref_lp_id]

    # decide value of minimum consumption package
    min_cons_val = hh.N_B_min * pref_bp.p[end]

    if min_cons_val > hh.I[end] + hh.S[end]
        B = hh.I[end] + hh.S[end]
        push!(hh.B, B)
    else
        B = min((1-hh.s) * hh.I[end], hh.I[end] + hh.S[end])
        push!(hh.B, B)
    end

    # println(hh.B[end])

    # TODO utility determination still has to happen

    # determine prices including sales taxes
    p_bp_τˢ = pref_bp.p[end] * (1 + τˢ)
    p_lp_τˢ = pref_lp.p[end] * (1 + τˢ)

    if min_cons_val > hh.I[end] + hh.S[end]

        U_B = rand(Uniform(0, 1))
        U_L = rand(Uniform(U_B, 1))

        α = U_B / (U_B + U_L)
        β = 1 - α

        N_B = hh.N_B_min + (α/p_bp_τˢ) * (hh.B[end] - p_bp_τˢ * hh.N_B_min)
        N_L = (β/p_lp_τˢ) * (hh.B[end] - p_bp_τˢ * hh.N_B_min)
        
    else
        N_B = (hh.I[end] + hh.S[end]) / p_bp_τˢ
        N_L = 0
    end

    return N_B, N_L
end


function update_sat_req_wage_hh!(
    hh::Household, 
    ϵ::Float64, 
    UB :: Float64
    )

    # Update satisfying wage as wage level over 4 periods
    # TODO: figure out if this should be wage or income
    if length(hh.w) > 4
        hh.wˢ = mean(hh.w[end-4:end])
    else
        hh.wˢ = mean(hh.w)
    end

    if hh.employed
        hh.wʳ = hh.wʳ * (1 + ϵ)
    else
        hh.wʳ = max(UB/hh.L, hh.wˢ)
    end

end


"""
Lets households get income, either from UB or wage
"""
function get_income_hh!(
    hh::Household, 
    amount :: Float64
    )

    push!(hh.I, amount)
    hh.C += amount
end


function get_fired_hh!(hh::Household)
    hh.employed = false
    hh.employer = 0
end


"""
Lets employee be hired, saves employer id and new earned wage.
"""
function get_hired_hh!(
    hh::Household, 
    wᴼ::Float64,
    employer::Int,
    )

    hh.employed = true
    hh.employer = employer
    hh.T_unemp = 0
    push!(hh.w, wᴼ)
end