mutable struct Household <: AbstractAgent
    id :: Int                   # global id
    hh_id :: Int                # hh id
    employed :: Bool            # is employed
    employer                    # employer
    I :: Array{Float64}         # hist income
    Iᵉ :: Float64               # expected income
    L :: Float64                # labor units in household
    S :: Array{Float64}         # total savings
    Sᵈ :: Float64               # desired savings
    s :: Float64                # savings rate
    B :: Array{Float64}         # budget
    C :: Float64                # cash
    N_B_min :: Float64          # substistence level of basic goods
    w :: Array{Float64}         # wage
    wˢ :: Float64               # satisfying wage
    wʳ :: Float64               # requested wage
    wᵉ :: Float64               # expected wage
    ωI :: Float64               # memory param income expectation
    bg                          # prefered basic good provider
    lg                          # prefered luxury good provider
    T_unemp :: Int              # time periods unemployed
end

function initialize_hh(id, hh_id, employed)
    hh = Household(
        id,                     # global id
        hh_id,                  # household id
        employed,               # bool: employed
        nothing,                # employer
        [100],                  # I: hist income
        100,                    # Iᵉ: exp income
        100,                    # L: labor units
        [100],                  # S: total savings
        0,                      # Sᵈ: desired savings
        0,                      # s: savings rate
        [0.0],                  # B: budget
        100,                    # C: cash
        1000,                   # N_B_min: substistence level of basic goods
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


function pick_cp_hh!(hh, all_bp, all_lp)
    
    # TODO: do this in a more sophisticated way
    hh.bg = sample(all_bp)
    hh.lg = sample(all_lp) 

end


function compute_exp_income_hh!(hh, P_HU, P_UU, UB)

    # determine income expectation for employed workers
    if hh.employed

        # update expected wage
        # TODO: find a way to reset this when unemployed
        if length(hh.w) > 1
            hh.wᵉ = hh.ωI * hh.wᵉ + (1-hh.ωI) * (2*hh.w[end] - hh.w[end-1])
        end

        P_UE = hh.employer.P_FE * (1 - P_HU)

        hh.Iᵉ = P_UE * UB + (1 - P_UE) * hh.wᵉ

    # determine income expectation for unemployed workers
    else

        hh.Iᵉ = P_UU * UB + (1 - P_UU) * hh.wʳ

    end

end


"""
Determines savings rate s
"""
function set_savingsrate_hh!(hh, avg_T_unemp, UB)

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
        hh.s = (hh.Sᵈ - hh.S[end]) / (hh.I[end] + 3*hh.Iᵉ)
    else

        hh.s = (B̄ - UB) / UB

    end

end


"""
Determines the optimal consumption package
    - Choose which goods to buy
    - Choose optimal consumption package for both goods
"""
function set_cons_package_hh!(hh, τˢ)

    # decide value of minimum consumption package
    min_cons_val = hh.N_B_min * hh.bg.p[end]

    if min_cons_val > hh.I[end] + hh.S[end]
        B = hh.I[end] + hh.S[end]
        push!(hh.B, B)
    else
        B = min((1-hh.s[end]) * hh.I[end], hh.I[end] + hh.S[end])
        push!(hh.B, B)
    end

    # TODO utility determination still has to happen
    p_bg_τˢ = hh.bg.p[end] * (1 - τˢ)
    p_lg_τˢ = hh.lg.p[end] * (1 - τˢ)

    U_B = rand(Uniform(0, 1))
    U_L = rand(Uniform(U_B, 1))

    α = U_B / (U_B + U_L)
    β = 1 - α

    N_B = hh.N_B_min + (α/p_bg_τˢ) * (hh.B[end] - p_bg_τˢ * hh.N_B_min)
    N_L = (β/p_lg_τˢ) * (hh.B[end] - p_bg_τˢ * hh.N_B_min)

    return N_B, N_L
end


function update_sat_req_wage_hh!(hh, ϵ :: Float64, UB :: Float64)

    # update satisfying wage as wage level over 4 periods
    if length(hh.w) > 4
        hh.wˢ = mean(hh.w[end-4:end])
    end

    if hh.employed
        hh.wʳ = hh.wʳ * (1 + ϵ)
    else
        hh.wʳ = max(UB/hh.L, hh.wˢ)
    end

end


function get_fired_hh!(hh)
    hh.employed = false
    hh.employer = nothing
end


function get_hired_hh!(hh, p)
    hh.employed = true
    hh.employer = p
    hh.T_unemp = 0
    push!(hh.w, p.wᴼ)
end