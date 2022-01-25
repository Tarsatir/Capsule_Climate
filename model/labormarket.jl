mutable struct LaborMarket
    employed :: Array{AbstractAgent}            # array of employed households
    unemployed :: Array{AbstractAgent}          # array of unemployed households
    E :: Float64                                # unemployment rate
    n_rounds :: Int                             # number of rounds in matching process
end

function initialize_labormarket()
    labormarket_struct = LaborMarket(
        [],
        [],
        0.05,
        1
    )
    return labormarket_struct
end

function update_unemploymentrate_lm(LM)
    LM.E = length(LM.unemployed) / (length(LM.employed) + length(LM.unemployed))
end

"""
Gives all producers a share of the labor force after being initialized
"""
function spread_employees_lm!(LM, all_cp, all_kp)

    i = 1
    for cp in all_cp
        employees = LM.employed[i:i+8]
        cp.Emp = employees
        cp.L = sum(map(hh -> hh.L, employees))
        i += 9
    end

    for kp in all_kp
        employees = LM.employed[i:i+9]
        kp.Emp = employees
        kp.L = sum(map(hh -> hh.L, employees))
        i += 9
    end

end

function matching_lm(labormarket_struct, all_cp, all_kp)

    # get all applicant workers
    # TODO: let employed workers also apply for jobs
    Lᵃ = labormarket_struct.unemployed

    for n_round in labormarket_struct.n_rounds

        

    end
end