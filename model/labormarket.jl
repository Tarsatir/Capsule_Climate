mutable struct LaborMarket
    employed :: Array{AbstractAgent}            # array of employed households
    unemployed :: Array{AbstractAgent}          # array of unemployed households
    E :: Float64                                # unemployment rate
end

function initialize_labormarket()
    labormarket_struct = LaborMarket(
        [],
        [],
        0.05
    )
    return labormarket_struct
end

function update_unemploymentrate_lm(LM)
    LM.E = length(LM.unemployed) / (length(LM.employed) + length(LM.unemployed))
end

function spread_employees_lm!(LM, all_cp, all_kp)

end

function matching_lm(labormarket_struct, all_cp, all_kp)
    labor_supply = sum(map(x -> x.L, labormarket_struct.employed))
    println(labor_supply)
    labor_demand = 0
    labor_demand += sum(map(cp -> cp.ΔLᵈ, all_cp))
    labor_demand += sum(map(kp -> kp.ΔLᵈ, all_kp))

    for cp in all_cp
        println(cp.ΔLᵈ)
    end

    println("yeet")

    for kp in all_kp
        println(kp.ΔLᵈ)
    end
    
    println(labor_demand)
end