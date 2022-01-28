mutable struct LaborMarket
    employed :: Array{AbstractAgent}            # array of employed households
    unemployed :: Array{AbstractAgent}          # array of unemployed households
    E :: Float64                                # unemployment rate
    n_rounds :: Int                             # number of rounds in matching process
    avg_T_unemp :: Float64                      # average time periods of unemployment
    P_HU :: Float64                             # hired workers as fraction of total unemployed workers
    P_UU :: Float64                             # probability of staying unemployed
end

function initialize_labormarket()
    labormarket_struct = LaborMarket(
        [],
        [],
        0.1,
        3,
        0,
        0,
        0
    )
    return labormarket_struct
end


"""
Defines the whole labor market process

- Lets producers fire excess workers
- Matches unemployed workers with producers looking for workers
- Updates unemployment rate

"""
function labormarket_process!(labormarket_struct, all_cp, all_kp, ϵ :: Float64, UB :: Float64)

    # get sets of firing and hiring producers
    firing_producers = []
    hiring_producers = []

    update_unemploymentrate_lm(labormarket_struct)
    println(labormarket_struct.E)

    for p in vcat(all_cp, all_kp)
        if p.ΔLᵈ > 0
            push!(hiring_producers, p)
        elseif p.ΔLᵈ < -100
            push!(firing_producers, p)
        end
    end

    # let producers fire excess workers
    fire_workers!(labormarket_struct, firing_producers)

    # update wage parameters households
    for hh in vcat(labormarket_struct.employed, labormarket_struct.unemployed)
        update_sat_req_wage_hh!(hh, ϵ, UB)
    end

    update_unemploymentrate_lm(labormarket_struct)
    println(labormarket_struct.E)

    # labor market matching process
    matching_lm(labormarket_struct, all_cp, all_kp)

    # update the unemployment rate
    update_unemploymentrate_lm(labormarket_struct)
    println(labormarket_struct.E)

end


function update_unemploymentrate_lm(labormarket_struct)
    labormarket_struct.E = (length(labormarket_struct.unemployed) / (length(labormarket_struct.employed) + length(labormarket_struct.unemployed)))
end


"""
Gives all producers a share of the labor force after being initialized
"""
function spread_employees_lm!(labormarket_struct, all_cp, all_kp)

    i = 1
    for cp in all_cp
        employees = labormarket_struct.employed[i:i+8]
        for emp in employees
            emp.employer = cp
        end
        cp.Emp = employees
        cp.L = sum(map(hh -> hh.L, employees))
        i += 9
    end

    for kp in all_kp
        employees = labormarket_struct.employed[i:i+8]
        for emp in employees
            emp.employer = kp
        end
        kp.Emp = employees
        kp.L = sum(map(hh -> hh.L, employees))
        i += 9
    end

end

function fire_workers!(labormarket_struct, firing_producers)

    all_fired_workers = []

    # choose who gets fired
    for p in firing_producers
        append!(all_fired_workers, fire_excess_workers_p!(p))
    end

    # change employment status for households
    for hh in all_fired_workers
        get_fired_hh!(hh)
    end

    # update employed and unemployed lists
    filter!(e -> e ∉ all_fired_workers, labormarket_struct.employed)
    append!(labormarket_struct.unemployed, all_fired_workers)

end


function matching_lm(labormarket_struct, all_cp, all_kp)

    # get all applicant workers
    # TODO: let employed workers also apply for jobs
    Lᵃ = labormarket_struct.unemployed

    n_unemployed = length(labormarket_struct.unemployed)

    # get producers that want to hire workers
    all_p = []
    for p in vcat(all_cp, all_kp)
        if p.ΔLᵈ > 0
            push!(all_p, p)
        end
    end

    n_hired = 0

    for n_round in 1:labormarket_struct.n_rounds

        # loop over producers
        shuffle!(all_p)
        for p in all_p

            # TODO do this in a more sophisticated way
            n_sample = min(10, length(Lᵃ))

            # make queue of job-seeking households
            Lˢ = sample(Lᵃ, n_sample, replace=false)

            # only select households with a low enough reservation wage
            Lᵈ = sort(Lˢ, by = l -> l.wʳ)

            # hire workers until demand is met or no more workers available
            while p.ΔLᵈ > 0 && length(Lᵈ) > 0

                # hire worker
                l = Lᵈ[1]
                get_hired_hh!(l, p)
                hire_worker_p!(p, l)

                n_hired += 1

                # delete household from seeking lists
                filter!(w -> w != l, Lᵈ)
                filter!(w -> w != l, Lᵃ)

            end
        end
    end

    # increase unemployment time for households
    n_longtermunemp = 0
    for hh in labormarket_struct.unemployed
        hh.T_unemp += 1
        if hh.T_unemp >= 2
            n_longtermunemp += 1
        end
    end

    labormarket_struct.P_UU = n_longtermunemp / n_unemployed
    labormarket_struct.P_HU = n_hired / n_unemployed

end


function update_avg_T_unemp_lm(labormarket_struct)
    labormarket_struct.avg_T_unemp = mean(map(hh -> hh.T_unemp, labormarket_struct.unemployed))
end