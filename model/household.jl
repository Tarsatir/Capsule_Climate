mutable struct Household <: AbstractAgent
    id :: Int                   # global id
    hh_id :: Int                # hh id
    employed :: Bool            # is employed
    employer                    # employer
    I :: Array{Float64}         # hist income
    Iᵉ :: Float64               # expected income
    L :: Float64                # labor units in household
    S :: Array{Float64}         # total savings
    B :: Float64                # budget
    w :: Array{Float64}         # wage
    wˢ :: Float64               # satisfying wage
    wʳ :: Float64               # requested wage
    ωI :: Float64               # memory param income expectation
    bg                          # prefered basic good provider
    lg                          # prefered luxury good provider
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
        0.0,                    # B: budget
        [1.0],                  # w: wage
        1.0,                    # wˢ: satisfying wage
        1.0,                    # wʳ: requested wage
        0.5,                    # ωI: memory param income exp
        nothing,                # bg: prefered basic good provider
        nothing                 # lg: prefered luxury good provider
    )
    return hh
end


function pick_cp_hh(hh, all_bp, all_lp)
    


end


"""
Determines consumption budget B
    - Determine income at time t
    - Estimate future income
    - Determine savings rate
    - Set consumption budget
"""
function set_budget_hh!(hh, UB, U, r)
    
    # determine income
    # Iₜ = UB
    # if (hh.employed)
    #     Iₜ = hh.w * hh.L
    # end
    # push!(hh.I, Iₜ)

    # hh.Iᵉ = compute_exp_income_hh(hh, U, r)

    # # determine savings rate
    # s = (Iₜ - Iᵉ) / Iₜ

    # TODO: better way to set up savings rate
    s = 0.1

    # Bₜ = 

end


"""
Determines the optimal consumption package
    - Choose which goods to buy
    - Choose optimal consumption package for both goods
"""
function set_cons_package_hh!(hh)

end


function compute_exp_income_hh(hh, U, r)

    ξ = 0
    if (length(U) > 2)
        if U[end] > U[end-1]
            ξ = rand(0, U[end] - U[end-1])
        else
            ξ = rand(U[end] - U[end-1], 0)
        end
    end

    Iᵉ = ωI * hh.Iᵉ + (1 - ωI) * (2 * hh.I[end] - hh.I[end-1]) + ξ * hh.I[end] + r[end] * hh.S[end]
    return Iᵉ
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
    push!(hh.w, p.wᴼ)
end