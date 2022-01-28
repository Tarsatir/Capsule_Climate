
mutable struct MacroEconomy
    GDP :: Array{Float64}       # GDP over time
    CPI :: Array{Float64}       # inflation over time
    C :: Array{Float64}         # aggregate consumption over time
    w :: Array{Float64}         # wage over time
    U :: Array{Float64}         # unemployment over time
    AB :: Array{Float64}        # average labor productivity over time
    l :: Array{Float64}         # unfilled demand over time
    E_bar :: Array{Float64}     # average competetiveness over time
    r :: Array{Float64}         # interest rate over time
    Ls :: Array{Float64}        # labor supply over time
    Ld :: Array{Float64}        # labor demand over time
end


function initialize_macro()
    macro_struct = MacroEconomy(
        [],                     # GDP over time
        [],                     # CPI: inflation rate
        [],                     # aggregate consumption
        [],                     # wage
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

function update_macro_stats(macro_struct, all_hh, all_cp, all_kp)

    # compute GDP and add to array
    total_I = sum(map(hh -> hh.I[end], all_hh))
    total_Π = sum(map(cp -> cp.Π[end], all_cp))
    total_Π += sum(map(kp -> kp.Π[end], all_kp))

    GDP = total_I + total_Π
    push!(macro_struct.GDP, GDP)

end


function update_labor_stats(macro_struct, labormarket_struct)

end
