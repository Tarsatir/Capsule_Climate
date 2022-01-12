
mutable struct MacroEconomy
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
        [],                     # CPI: inflation rate
        [],                     # aggregate consumption
        [2],                    # wage
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


function compute_average_competitiveness!(macro_struct, consumer_good_producers)
    E_bar_t = sum(map(x -> x.E[end] * x.f[end], consumer_good_producers))
    push!(E_bar, macro_struct.E_bar_t)

    return E_bar_t
end


function update_wage!(macro_struct, global_param)
    ΔAB = AB[end] - AB[end-1]
    ΔCPI = CPI[end] - CPI[end-1]
    ΔU = U[end] - U[end-1]
    w_t = w[end] + (1 + global_param.ψ1*(ΔAB/AB[end-1]) 
                      + global_param.ψ2*(ΔCPI/CPI[end-1]) 
                      + global_param.ψ3*(ΔU/U[end-1]))
    push!(w_t, macro_struct.w)
end
