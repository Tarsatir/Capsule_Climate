mutable struct Household <: AbstractAgent

    id :: Int                   # global id

    # Employment variables
    employed :: Bool            # is employed
    employer_id :: Int          # id of employer
    L :: Float64                # labor units in household
    w :: Vector{Float64}        # wage
    wˢ :: Float64               # satisfying wage
    wʳ :: Float64               # requested wage
    T_unemp :: Int              # time periods unemployed

    # Income and wealth variables
    I :: Vector{Float64}        # hist income
    Iᵀ :: Vector{Float64}       # hist taxed income
    S :: Vector{Float64}        # total savings
    s :: Float64                # savings rate
    W :: Vector{Float64}        # wealth or cash on hand
    Wʳ :: Vector{Float64}       # real wealth or cash on hand

    # Consumption variables
    C :: Vector{Float64}        # budget
    N_B_min :: Float64          # substistence level of basic goods
    bp :: Vector{Int}           # connected cp basic goods
    lp :: Vector{Int}           # connected cp luxury goods
    P̄ :: Float64                # weighted average price of bp
    c_L :: Float64              # share of income used to buy luxury goods

end

function initialize_hh(
    id::Int,
    )::Household

    hh = Household(
        id,                     # global id

        false,                  # bool: employed
        0,                      # id of employer
        100,                    # L: labor units
        [1.0],                  # w: wage
        1.0,                    # wˢ: satisfying wage
        1.0,                    # wʳ: requested wage
        0,                      # T_unemp: time periods unemployed

        [],                     # I: hist income
        [],                     # Iᵀ: hist taxed income
        [10],                   # S: total savings
        0,                      # s: savings rate
        [],                     # W: wealth or cash on hand
        [],                     # Wʳ: real wealth or cash on hand

        [0.0],                  # C: budget
        10,                     # N_B_min: substistence level of basic goods
        Vector{Int}(),          # all_cp_B: connected cp basic goods
        Vector{Int}(),          # all_cp_L: connected cp luxury goods
        0,                      # P̄: weighted average price of bp
        0.5,                    # c_L: share of budget used to buy luxury goods
    )
    return hh
end


"""
Uniformly samples basic and luxury good producers to be in trading network.
"""
function select_bp_lp_hh!(
    hh::Household,
    all_bp::Vector{Int},
    all_lp::Vector{Int},
    n_bp::Int,
    n_lp::Int
    )

    selected_bp = sample(all_bp, n_bp)
    hh.bp = selected_bp

    selected_lp = sample(all_lp, n_lp)
    hh.lp = selected_lp
end


"""
Updates current wealth level W, equal to the current cash on hand
"""
function update_W_hh!(
    hh::Household
    )
    
    W = hh.Iᵀ[end] + hh.S[end]
    push!(hh.W, W)
end


"""
Computes average price level of available bp and lp
"""
function update_average_price_hh!(
    hh::Household,
    model::ABM
    )

    P̄_bp = mean(bp_id -> model[bp_id].p[end], hh.bp)
    P̄_lp = mean(lp_id -> model[lp_id].p[end], hh.lp)
    hh.P̄ = hh.c_L * P̄_lp + (1 - hh.c_L) * P̄_bp
end


"""
Updates real wealth level based on price division of last period
"""
function update_real_wealth_hh!(
    hh::Household
    )

    Wʳ = hh.W[end] / hh.P̄[end]
    push!(hh.Wʳ, Wʳ)
end


"""
Defines logistic function that determines share of each good type
"""
function update_share_goodtypes_hh!(
    hh::Household,
    c_L_max::Float64,
    a_σ::Float64,
    b_σ::Float64
    )

    hh.c_L = c_L_max / (1 + exp(-(hh.Wʳ[end]/a_σ - b_σ)))
end


"""
Computes consumption budget
"""
function compute_consumption_budget_hh!(
    hh::Household,
    α_cp::Float64
    )

    C = min((hh.Wʳ[end])^α_cp, hh.Wʳ[end])
    push!(hh.C, C)
end


"""
Sets consumption budget based on current wealth level
"""
function set_consumption_budget_hh!(
    hh::Household,
    a_σ::Float64,
    b_σ::Float64,
    α_cp::Float64,
    model::ABM
    )

    # Update average price level of bp and lp
    update_average_price_hh!(hh, model)

    # Update real wealth level
    update_real_wealth_hh!(hh)

    # Update share of goods to bg and lg
    update_share_goodtypes_hh!(hh, c_L_max, a_σ, b_σ)

    # Compute consumption budget
    compute_consumption_budget_hh!(hh, α_cp)
end

# function pick_cp_hh!(
#     hh::Household, 
#     supplying_bp::Vector{Int}, 
#     supplying_lp::Vector{Int}
#     )
    
#     # TODO: do this in a more sophisticated way
#     hh.pref_bp_id = sample(supplying_bp)
#     hh.pref_lp_id = sample(supplying_lp)

# end


# """
# Computes the expected income based on perceived probabilities.
# """
# function compute_exp_income_hh!(
#     hh::Household, 
#     P_HU::Float64, 
#     P_UU::Float64, 
#     UB::Float64,
#     model::ABM
#     )

#     # determine income expectation for employed workers
#     if hh.employed

#         # update expected wage
#         # TODO: find a way to reset this when unemployed
#         if length(hh.w) > 1
#             hh.wᵉ = hh.ωI * hh.wᵉ + (1-hh.ωI) * (2*hh.w[end] - hh.w[end-1])
#         end

#         P_UE = model[hh.employer_id].P_FE * (1 - P_HU)

#         hh.Iᵉ = P_UE * UB + (1 - P_UE) * hh.wᵉ

#     # determine income expectation for unemployed workers
#     else

#         hh.Iᵉ = P_UU * UB + (1 - P_UU) * hh.wʳ

#     end
# end


# """
# Determines savings rate s
# """
# function set_savingsrate_hh!(
#     hh::Household, 
#     avg_T_unemp::Float64, 
#     UB::Float64
#     )

#     # determine average budget over last year or over available information
#     if length(hh.B[end]) >= 5
#         B̄ = mean(hh.B[end-4, end])
#     else
#         B̄ = mean(hh.B)
#     end
    
#     if hh.employed
#         # determine desired level of savings
#         hh.Sᵈ = avg_T_unemp * (B̄ - UB)

#         # determine savings rate
#         s = max((hh.Sᵈ - hh.S[end]) / (hh.I[end] + 3*hh.Iᵉ), -hh.S[end]/hh.I[end])
#     else

#         s = (B̄ - UB) / UB

#     end

#     # TODO: find a solution for this
#     if isnan(s)
#         hh.s = 0
#     else
#         hh.s = s
#     end

#     # println(hh.s, " ", hh.I[end], " ", 3*hh.Iᵉ)

# end


# """
# Determines the optimal consumption package
#     - Choose which goods to buy
#     - Choose optimal consumption package for both goods
# """
# function set_cons_package_hh!(
#     hh::AbstractAgent, 
#     τˢ::Float64,
#     model::ABM
#     )::Tuple{Float64, Float64}

#     pref_bp = model[hh.pref_bp_id]
#     pref_lp = model[hh.pref_lp_id]

#     # decide value of minimum consumption package
#     min_cons_val = hh.N_B_min * pref_bp.p[end]

#     if min_cons_val > hh.I[end] + hh.S[end]
#         B = hh.I[end] + hh.S[end]
#         push!(hh.B, B)
#     else
#         B = min((1-hh.s) * hh.I[end], hh.I[end] + hh.S[end])
#         push!(hh.B, B)
#     end

#     # println(hh.B[end])

#     # TODO utility determination still has to happen

#     # determine prices including sales taxes
#     p_bp_τˢ = pref_bp.p[end] * (1 + τˢ)
#     p_lp_τˢ = pref_lp.p[end] * (1 + τˢ)

#     if min_cons_val > hh.I[end] + hh.S[end]

#         U_B = rand(Uniform(0, 1))
#         U_L = rand(Uniform(U_B, 1))

#         α = U_B / (U_B + U_L)
#         β = 1 - α

#         N_B = hh.N_B_min + (α/p_bp_τˢ) * (hh.B[end] - p_bp_τˢ * hh.N_B_min)
#         N_L = (β/p_lp_τˢ) * (hh.B[end] - p_bp_τˢ * hh.N_B_min)
        
#     else
#         N_B = (hh.I[end] + hh.S[end]) / p_bp_τˢ
#         N_L = 0
#     end

#     return N_B, N_L
# end


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


function set_unemployed_hh!(
    hh::Household
    )

    hh.employed = false
    hh.employer_id = 0
end


"""
Lets employee be hired, saves employer id and new earned wage.
"""
function set_employed_hh!(
    hh::Household, 
    wᴼ::Float64,
    employer_id::Int,
    )

    hh.employed = true
    hh.employer_id = employer_id
    hh.T_unemp = 0
    push!(hh.w, wᴼ)
end