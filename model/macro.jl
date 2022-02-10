
mutable struct MacroEconomy
    GDP :: Array{Float64}       # GDP over time
    GDP_I :: Array{Float64}     # income share of GDP over time
    GDP_Π_cp :: Array{Float64}  # profit share of GDP of cp over time
    GDP_Π_kp :: Array{Float64}  # profit share of GDP of kp over time
    CPI :: Array{Float64}       # inflation over time
    C :: Array{Float64}         # aggregate consumption over time
    w̄_avg :: Array{Float64}     # average wage over time
    w̄_std :: Array{Float64}     # std of wage over time
    Ī_avg :: Array{Float64}     # average income over time
    Ī_std :: Array{Float64}     # std of income over time
    B̄_avg :: Array{Float64}     # average of budget over time
    B̄_std :: Array{Float64}     # std of budget over time
    U :: Array{Float64}         # unemployment over time
    Exp_UB :: Array{Float64}    # total spending on UB
    AB :: Array{Float64}        # average labor productivity over time
    l :: Array{Float64}         # unfilled demand over time
    E_bar :: Array{Float64}     # average competetiveness over time
    r :: Array{Float64}         # interest rate over time
    Ls :: Array{Float64}        # labor supply over time
    Ld :: Array{Float64}        # labor demand over time

    # Changes in savings rate
    s̄_avg :: Array{Float64}     # average savings rate over time
    s̄_std :: Array{Float64}     # std of savings over time

    # Changes in labor demand
    ΔL̄_avg :: Array{Float64}    # average desired labor change
    ΔL̄_std :: Array{Float64}    # std desired labor change
    ΔL̄_cp_avg :: Array{Float64} # average desired over time for cp
    ΔL̄_cp_std :: Array{Float64} # std desired labor change for cp
    ΔL̄_kp_avg :: Array{Float64} # ' ' for kp
    ΔL̄_kp_std :: Array{Float64}
end


function initialize_macro()
    macro_struct = MacroEconomy(
        [],                     # GDP over time
        [],                     # income share of GDP
        [],                     # profit share of GDP of cp over time
        [],                     # profit share of GDP of kp over time
        [],                     # CPI: inflation rate
        [],                     # aggregate consumption
        [],                     # average of wage over time
        [],                     # std of wage over time
        [],                     # average income over time
        [],                     # std of income over time
        [],                     # average budget over time
        [],                     # std of budget over time
        [],                     
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        []
    )
    return macro_struct
end


"""
Updates macro stats after each time step
"""
function update_macro_timeseries(
    macro_struct, 
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    all_kp::Vector{Int}, 
    E::Float64, 
    Exp_UB::Float64,
    model::ABM
    )

    # Retrieve agent structs
    all_hh_str = map(hh_id -> model[hh_id], all_hh)
    all_cp_str = map(cp_id -> model[cp_id], all_cp)
    all_kp_str = map(kp_id -> model[kp_id], all_kp)

    # Compute GDP and add to array
    compute_GDP!(all_hh_str, all_cp_str, all_kp_str, macro_struct)

    ΔL̄_avg = mean(map(p -> p.ΔLᵈ, vcat(all_cp_str, all_kp_str)))
    ΔL̄_std = std(map(p -> p.ΔLᵈ, vcat(all_cp_str, all_kp_str)))
    push!(macro_struct.ΔL̄_avg, ΔL̄_avg)
    push!(macro_struct.ΔL̄_std, ΔL̄_std)

    ΔL̄_cp_avg = mean(map(cp -> cp.ΔLᵈ, all_cp_str))
    ΔL̄_cp_std = std(map(cp -> cp.ΔLᵈ, all_cp_str))
    push!(macro_struct.ΔL̄_cp_avg, ΔL̄_cp_avg)
    push!(macro_struct.ΔL̄_cp_std, ΔL̄_cp_std)

    ΔL̄_kp_avg = mean(map(kp -> kp.ΔLᵈ, all_kp_str))
    ΔL̄_kp_std = std(map(kp -> kp.ΔLᵈ, all_kp_str))
    push!(macro_struct.ΔL̄_kp_avg, ΔL̄_kp_avg)
    push!(macro_struct.ΔL̄_kp_std, ΔL̄_kp_std)

    push!(macro_struct.U, E)

    push!(macro_struct.Exp_UB, Exp_UB)

    # compute average savings rate
    s̄_avg = mean(map(hh -> hh.s, all_hh_str))
    s̄_std = std(map(hh -> hh.s, all_hh_str))
    # println(map(hh -> hh.s, all_hh_str))
    push!(macro_struct.s̄_avg, s̄_avg)
    push!(macro_struct.s̄_std, s̄_std)

    # TODO: also include kp
    w̄_avg = mean(map(cp -> cp.w[end], all_cp_str))
    w̄_std = std(map(cp -> cp.w[end], all_cp_str))
    push!(macro_struct.w̄_avg, w̄_avg)
    push!(macro_struct.w̄_std, w̄_std)

    Ī_avg = mean(map(hh -> hh.I[end], all_hh_str))
    Ī_std = std(map(hh -> hh.I[end], all_hh_str))
    push!(macro_struct.Ī_avg, Ī_avg)
    push!(macro_struct.Ī_std, Ī_std)

    D̄ = mean(map(cp -> cp.D[end], all_cp_str))
    println(D̄)
end


"""
Computes GDP based on income of separate sectors, computes partial incomes of sectors
"""
function compute_GDP!(
    all_hh_str::Vector{Household},
    all_cp_str::Vector{ConsumerGoodProducer},
    all_kp_str::Vector{CapitalGoodProducer},
    macro_struct
    )
    total_I = sum(map(hh -> hh.I[end], all_hh_str))
    push!(macro_struct.GDP_I, total_I)

    total_Π_cp = sum(map(cp -> cp.Π[end], all_cp_str))
    push!(macro_struct.GDP_Π_cp, total_Π_cp)
    
    total_Π_kp = sum(map(kp -> kp.Π[end], all_kp_str))
    push!(macro_struct.GDP_Π_kp, total_Π_kp)

    GDP = total_I + total_Π_cp + total_Π_kp
    push!(macro_struct.GDP, GDP)
end

function update_labor_stats(macro_struct, labormarket_struct)

end
