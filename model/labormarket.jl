mutable struct LaborMarket
    employed :: Array{AbstractAgent}            # array of employed households
    unemployed :: Array{AbstractAgent}          # array of unemployed households
    # employed :: Array{Int}
    # unemployed :: Array{Int}
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
        1,
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
function labormarket_process!(labormarket_struct, all_cp :: Array{AbstractAgent}, all_kp :: Array{AbstractAgent}, ϵ :: Float64, UB :: Float64)

    # get sets of firing and hiring producers
    firing_producers = []
    hiring_producers = []

    update_unemploymentrate_lm(labormarket_struct)
    println("E 1: ", labormarket_struct.E[end])

    for p in vcat(all_cp, all_kp)
        if p.ΔLᵈ > 0
            push!(hiring_producers, p)
        elseif p.ΔLᵈ < -100
            push!(firing_producers, p)
        end
    end

    # let producers fire excess workers
    println("f1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    fire_workers!(labormarket_struct, firing_producers)
    println("f2 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

    # update wage parameters households
    for hh in vcat(labormarket_struct.employed, labormarket_struct.unemployed)
        update_sat_req_wage_hh!(hh, ϵ, UB)
    end

    update_unemploymentrate_lm(labormarket_struct)
    println("E 2: ", labormarket_struct.E)

    # labor market matching process
    println("m1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    matching_lm(labormarket_struct, hiring_producers)
    println("m1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    
    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

    # update the unemployment rate
    update_unemploymentrate_lm(labormarket_struct)
    println("E 3: ", labormarket_struct.E)

end


function update_unemploymentrate_lm(labormarket_struct)
    labormarket_struct.E = (length(labormarket_struct.unemployed) / (length(labormarket_struct.employed) + length(labormarket_struct.unemployed)))
end


"""
Gives all producers a share of the labor force after being initialized
"""
function spread_employees_lm!(labormarket_struct, all_hh, all_cp, all_kp, n_init_emp_cp, n_init_emp_kp)

    # n_init_emp_cp = 11
    # n_init_emp_kp = 3

    i = 1
    for cp in all_cp
        employees = all_hh[i:i+n_init_emp_cp-1]
        for emp in employees
            # emp.employer = cp
            push!(labormarket_struct.employed, emp)
            get_hired_hh!(emp, cp)
            hire_worker_p!(cp, emp)
            
        end
        cp.Emp = employees
        cp.L = sum(map(hh -> hh.L, employees))
        i += n_init_emp_cp
    end

    for kp in all_kp
        employees = all_hh[i:i+n_init_emp_kp-1]
        for emp in employees
            # emp.employer = kp
            push!(labormarket_struct.employed, emp)
            get_hired_hh!(emp, kp)
            hire_worker_p!(kp, emp)
        end
        kp.Emp = employees
        kp.L = sum(map(hh -> hh.L, employees))
        i += n_init_emp_kp
    end

    # other households are pushed into unemployment
    for hh in all_hh[i:end]
        get_fired_hh!(hh)
        push!(labormarket_struct.unemployed, hh)
    end

end


function fire_workers!(labormarket_struct, firing_producers)

    all_fired_workers = []

    # choose who gets fired
    for p in firing_producers
        append!(all_fired_workers, fire_excess_workers_p!(p))
    end

    # update employed and unemployed lists
    filter!(e -> e ∉ all_fired_workers, labormarket_struct.employed)

    # change employment status for households
    for hh in all_fired_workers
        get_fired_hh!(hh)
    end

    append!(labormarket_struct.unemployed, all_fired_workers)

end


function matching_lm(labormarket_struct, hiring_producers)

    # get all applicant workers
    # TODO: let employed workers also apply for jobs
    # Lᵃ = labormarket_struct.unemployed

    n_unemployed = length(labormarket_struct.unemployed)

    n_hired = 0

    for n_round in 1:labormarket_struct.n_rounds

        # loop over producers
        shuffle!(hiring_producers)
        for p in hiring_producers

            if p.ΔLᵈ > 0

                # TODO do this in a more sophisticated way
                n_sample = min(10, length(labormarket_struct.unemployed))

                # make queue of job-seeking households
                Lˢ = sample(labormarket_struct.unemployed, n_sample, replace=false)

                # only select households with a low enough reservation wage
                Lᵈ = sort(Lˢ, by = l -> l.wʳ)

                to_be_hired = []
                w = p.w[end]
                for l in Lᵈ
                    # hire worker
                    l = Lᵈ[1]
                    # n_hired += 1

                    push!(to_be_hired, l)
                    w = l.wʳ
                    # println(l.wʳ)

                    # delete household from seeking lists
                    # filter!(hh -> hh ≠ l, Lᵈ)
                    # filter!(hh -> hh ≠ l, labormarket_struct.unemployed)
                    # push!(labormarket_struct.employed, l)
                end

                # add wage at which hired to wage list
                p.wᴼ = w
                # println(p.wᴼ)
                # push!(p.w, w)

                # hire workers until demand is met or no more workers available
                # while p.ΔLᵈ > 0 && length(Lᵈ) > 0
                for l in to_be_hired
                    
                    # delete household from seeking lists
                    filter!(hh -> hh ≠ l, Lᵈ)

                    filter!(hh -> hh ≠ l, labormarket_struct.unemployed)

                    get_hired_hh!(l, p)
                    hire_worker_p!(p, l)

                    push!(labormarket_struct.employed, l)

                    # hire worker
                    # l = Lᵈ[1]
                    

                    n_hired += 1

                end

                update_wage_level_p!(p)

            end

        end
    end

    # increase unemployment time for households
    # TODO put this in households
    n_longtermunemp = 0
    for hh in labormarket_struct.unemployed
        hh.T_unemp += 1
        if hh.T_unemp >= 2
            n_longtermunemp += 1
        end
    end

    labormarket_struct.P_UU = n_longtermunemp / n_unemployed
    # labormarket_struct.P_HU = n_hired / n_unemployed
    labormarket_struct.P_HU = 1 - labormarket_struct.P_UU

    # println(labormarket_struct.P_UU, " ", labormarket_struct.P_HU)

end


function update_avg_T_unemp_lm(labormarket_struct)
    if length(labormarket_struct.unemployed) > 0
        labormarket_struct.avg_T_unemp = mean(map(hh -> hh.T_unemp, labormarket_struct.unemployed))
    else
        labormarket_struct.avg_T_unemp = 0
    end
end