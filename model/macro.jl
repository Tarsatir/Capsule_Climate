
mutable struct MacroEconomy
    GDP :: Array{Float64}       # GDP over time
    GDP_I :: Array{Float64}     # income share of GDP over time
    GDP_Π_cp :: Array{Float64}  # profit share of GDP of cp over time
    GDP_Π_kp :: Array{Float64}  # profit share of GDP of kp over time
    CPI :: Array{Float64}       # inflation over time
    C :: Array{Float64}         # aggregate consumption over time
    w̄_avg :: Array{Float64}     # average wage over time
    w̄_std :: Array{Float64}     # std of wage over time
    U :: Array{Float64}         # unemployment over time
    AB :: Array{Float64}        # average labor productivity over time
    l :: Array{Float64}         # unfilled demand over time
    E_bar :: Array{Float64}     # average competetiveness over time
    r :: Array{Float64}         # interest rate over time
    Ls :: Array{Float64}        # labor supply over time
    Ld :: Array{Float64}        # labor demand over time
    s̄_avg :: Array{Float64}     # average savings rate over time
    s̄_std :: Array{Float64}     # std of savings over time
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

function update_macro_stats(macro_struct, all_hh, all_cp, all_kp, E)

    # compute GDP and add to array
    total_I = sum(map(hh -> hh.I[end], all_hh))
    push!(macro_struct.GDP_I, total_I)

    total_Π_cp = sum(map(cp -> cp.Π[end], all_cp))
    push!(macro_struct.GDP_Π_cp, total_Π_cp)
    
    total_Π_kp = sum(map(kp -> kp.Π[end], all_kp))
    push!(macro_struct.GDP_Π_kp, total_Π_kp)

    GDP = total_I + total_Π_cp + total_Π_kp
    push!(macro_struct.GDP, GDP)

    push!(macro_struct.U, E)

    # compute average savings rate
    s̄_avg = mean(map(hh -> hh.s, all_hh))
    s̄_std = std(map(hh -> hh.s, all_hh))
    # println(map(hh -> hh.s, all_hh))
    push!(macro_struct.s̄_avg, s̄_avg)
    push!(macro_struct.s̄_std, s̄_std)

    # TODO: also include kp
    w̄_avg = mean(map(cp -> cp.w[end], all_cp))
    w̄_std = std(map(cp -> cp.w[end], all_cp))
    push!(macro_struct.w̄_avg, w̄_avg)
    push!(macro_struct.w̄_std, w̄_std)
    

end


function update_labor_stats(macro_struct, labormarket_struct)

end
