
mutable struct MacroEconomy
    CPI :: Array{Float64}
    C :: Array{Float64}
    w :: Array{Float64}
    U :: Array{Float64}
    AB :: Array{Float64}
end

function update_wage(macro_struct, global_param)
    ΔAB = AB[end] - AB[end-1]
    ΔCPI = CPI[end] - CPI[end-1]
    ΔU = U[end] - U[end-1]
    w_t = w[end] + (1 + global_param.ψ1*(ΔAB/AB[end-1]) 
                      + global_param.ψ2*(ΔCPI/CPI[end-1]) 
                      + global_param.ψ3*(ΔU/U[end-1]))
end
