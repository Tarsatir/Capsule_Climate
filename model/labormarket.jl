mutable struct LaborMarket
    employed :: Vector{Int}           # array of employed households
    unemployed :: Vector{Int}          # array of unemployed households
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
Gives all producers a share of the labor force after being initialized.
Households that are not assigned to a producer are declared unemployed.
"""
function spread_employees_lm!(
    labormarket_struct, 
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    all_kp::Vector{Int}, 
    n_init_emp_cp::Int, 
    n_init_emp_kp::Int,
    model::ABM
    )

    i = 1
    for cp_id in all_cp
        cp = model[cp_id]
        employees = all_hh[i:i+n_init_emp_cp-1]
        for hh_id in employees
            hh = model[hh_id]
            get_hired_hh!(hh, 1.0, cp_id)
            hire_worker_p!(cp, hh)
            push!(labormarket_struct.employed, hh.id)
        end

        cp.Emp = employees
        cp.L = sum(map(hh_id -> model[hh_id].L, employees))
        i += n_init_emp_cp

    end

    for kp_id in all_kp
        kp = model[kp_id]
        employees = all_hh[i:i+n_init_emp_kp-1]
        for hh_id in employees
            hh = model[hh_id]
            get_hired_hh!(hh, 1.0, kp_id)
            hire_worker_p!(kp, hh)
            push!(labormarket_struct.employed, hh.id)
        end
        kp.Emp = employees
        kp.L = sum(map(hh_id -> model[hh_id].L, employees))
        i += n_init_emp_kp
    end

    # other households are pushed into unemployment
    for hh_id in all_hh[i:end]
        hh = model[hh_id]
        get_fired_hh!(hh)
        push!(labormarket_struct.unemployed, hh.id)
    end

end


"""
Defines the whole labor market process

- Lets producers fire excess workers
- Matches unemployed workers with producers looking for workers
- Updates unemployment rate

"""
function labormarket_process!(
    labormarket_struct,
    all_hh::Vector{Int}, 
    all_p::Vector{Int}, 
    ϵ :: Float64, 
    UB :: Float64,
    model::ABM
    )

    # get sets of firing and hiring producers
    firing_producers = Vector{Int}()
    hiring_producers = Vector{Int}()

    update_unemploymentrate_lm(labormarket_struct)
    println("E 1: ", labormarket_struct.E[end])

    for p_id in all_p
        if model[p_id].ΔLᵈ > 0
            push!(hiring_producers, p_id)
        elseif model[p_id].ΔLᵈ < -100
            push!(firing_producers, p_id)
        end
    end

    # let producers fire excess workers
    # println("f1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))
    fire_workers!(labormarket_struct, firing_producers, model)
    # println("f2 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))

    update_employment_lm(labormarket_struct, all_hh, all_p, model)

    # update wage parameters households
    for hh_id in all_hh
        update_sat_req_wage_hh!(model[hh_id], ϵ, UB)
    end

    update_unemploymentrate_lm(labormarket_struct)
    println("E 2: ", labormarket_struct.E)

    # labor market matching process
    # println("m1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))
    matching_lm(labormarket_struct, hiring_producers, all_p, all_hh, model)
    # println("m1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))


    
    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

    # update the unemployment rate
    update_unemploymentrate_lm(labormarket_struct)
    println("E 3: ", labormarket_struct.E)

end


function update_unemploymentrate_lm(labormarket_struct)
    labormarket_struct.E = (length(labormarket_struct.unemployed) / (length(labormarket_struct.employed) + length(labormarket_struct.unemployed)))
end


function fire_workers!(
    labormarket_struct, 
    firing_producers::Vector{Int}, 
    model::ABM
    )

    all_fired_workers = Vector{Int}()

    # choose who gets fired
    for p_id in firing_producers
        fired_workers = fire_excess_workers_p!(model[p_id], model)
        append!(all_fired_workers, fired_workers)
    end

    # println("yeet ", length(all_fired_workers))

    # update employed and unemployed lists
    filter!(hh_id -> hh_id ∉ all_fired_workers, labormarket_struct.employed)
    append!(labormarket_struct.unemployed, all_fired_workers)

    # println("unemp", labormarket_struct.unemployed)

    # change employment status for households
    for hh_id in all_fired_workers
        get_fired_hh!(model[hh_id])
    end

end


function matching_lm(
    labormarket_struct, 
    hiring_producers::Vector{Int},
    all_p::Vector{Int},
    all_hh::Vector{Int},
    model::ABM
    )

    # get all applicant workers
    # TODO: let employed workers also apply for jobs
    # Lᵃ = labormarket_struct.unemployed

    n_unemployed = length(labormarket_struct.unemployed)

    n_hired = 0

    for n_round in 1:labormarket_struct.n_rounds

        # loop over producers
        shuffle!(hiring_producers)
        for p_id in hiring_producers

            # Stop process if no unemployed left
            if length(labormarket_struct.unemployed) == 0
                return
            end

            p = model[p_id]

            if p.ΔLᵈ > 0

                # TODO do this in a more sophisticated way
                n_sample = min(10, length(labormarket_struct.unemployed))

                # make queue of job-seeking households
                Lˢ = sample(labormarket_struct.unemployed, n_sample, replace=false)

                # only select households with a low enough reservation wage
                Lᵈ = sort(Lˢ, by = hh_id -> model[hh_id].wʳ)

                to_be_hired = []
                # w = p.w[end]

                ΔL = p.ΔLᵈ

                for hh_id in Lᵈ
                    # hire worker
                    # hh_id = Lᵈ[1]
                    # n_hired += 1
                    
                    ΔL -= model[hh_id].L
                    push!(to_be_hired, hh_id)

                    if ΔL < 0
                        break
                    end
                    # w = model[hh_id].wʳ
                    # println(l.wʳ)

                    # delete household from seeking lists
                    # filter!(hh -> hh ≠ l, Lᵈ)
                    # filter!(hh -> hh ≠ l, labormarket_struct.unemployed)
                    # push!(labormarket_struct.employed, l)
                end

                # add wage at which hired to wage list
                # p.wᴼ = w
                p.wᴼ = model[to_be_hired[end]].wʳ
                # println(p.wᴼ)
                # push!(p.w, w)

                # hire workers until demand is met or no more workers available
                # while p.ΔLᵈ > 0 && length(Lᵈ) > 0
                for hh_id in to_be_hired
                    
                    # delete household from seeking lists
                    # filter!(hh -> hh ≠ hh_id, Lᵈ)

                    get_hired_hh!(model[hh_id], p.wᴼ, p_id)
                    hire_worker_p!(model[p_id], model[hh_id])

                    # filter!(hh -> hh ≠ hh_id, labormarket_struct.unemployed)
                    # push!(labormarket_struct.employed, hh_id)                    

                    n_hired += 1

                end

                update_wage_level_p!(p, model)

            end

            # TODO find a way not to have to do this every time
            update_employment_lm(labormarket_struct, all_hh, all_p, model)

        end

    end

    # increase unemployment time for households
    # TODO put this in households
    n_longtermunemp = 0
    for hh_id in labormarket_struct.unemployed
        model[hh_id].T_unemp += 1
        if model[hh_id].T_unemp >= 2
            n_longtermunemp += 1
        end
    end

    labormarket_struct.P_UU = n_longtermunemp / n_unemployed
    # labormarket_struct.P_HU = n_hired / n_unemployed
    labormarket_struct.P_HU = 1 - labormarket_struct.P_UU

    # println(labormarket_struct.P_UU, " ", labormarket_struct.P_HU)

end


function update_avg_T_unemp_lm(
    labormarket_struct, 
    model::ABM
    )

    if length(labormarket_struct.unemployed) > 0
        avg_T_unemp = mean(map(hh_id -> model[hh_id].T_unemp, labormarket_struct.unemployed))
        labormarket_struct.avg_T_unemp = avg_T_unemp
    else
        labormarket_struct.avg_T_unemp = 0
    end
end


function update_employment_lm(
    labormarket_struct,
    all_hh::Vector{Int},
    all_p::Vector{Int},
    model::ABM
    )

    employed = Vector{Int}()

    for p_id in all_p
        for hh_id in model[p_id].Emp
            push!(employed, hh_id)
        end
    end

    unemployed = filter(hh_id -> hh_id ∉ employed, all_hh)

    labormarket_struct.employed = employed
    labormarket_struct.unemployed = unemployed

end