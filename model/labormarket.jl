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
    fire_workers_lm!(labormarket_struct, firing_producers, model)
    # println("f2 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))

    # Update wage parameters households
    for hh_id in all_hh
        update_sat_req_wage_hh!(model[hh_id], ϵ, UB)
    end

    update_unemploymentrate_lm(labormarket_struct)
    println("E 2: ", labormarket_struct.E)

    # labor market matching process
    # println("m1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))
    matching_lm(labormarket_struct, hiring_producers, model)
    # println("m2 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].Emp), all_p)))

    # Update probabilities of long-term unemployment
    update_probs_lm!(labormarket_struct, model)
    
    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

    # update the unemployment rate
    update_unemploymentrate_lm(labormarket_struct)
    println("E 3: ", labormarket_struct.E)

    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

end


function update_unemploymentrate_lm(labormarket_struct)
    labormarket_struct.E = (length(labormarket_struct.unemployed) / (length(labormarket_struct.employed) + length(labormarket_struct.unemployed)))
end


"""
Employers select workers to fire, labor market aggregates are changed and employees 
are fired.
"""
function fire_workers_lm!(
    labormarket_struct, 
    firing_producers::Vector{Int}, 
    model::ABM
    )

    all_fired_workers = Vector{Int}()

    # choose who gets fired
    for p_id in firing_producers
        fired_workers = fire_excess_workers_p!(model[p_id], model)
        # println(fired_workers)
        append!(all_fired_workers, fired_workers)
    end

    # println(all_fired_workers)

    # update employed and unemployed lists
    update_firedworker_lm!(labormarket_struct, all_fired_workers)

    # change employment status for households
    for hh_id in all_fired_workers
        get_fired_hh!(model[hh_id])
    end
end


function matching_lm(
    labormarket_struct, 
    hiring_producers::Vector{Int},
    model::ABM
    )

    # TODO: let employed workers also apply for jobs

    for n_round in 1:labormarket_struct.n_rounds

        # loop over hiring producers producers
        for p_id in hiring_producers

            # Stop process if no unemployed left
            if length(labormarket_struct.unemployed) == 0
                return
            end

            p = model[p_id]

            if p.ΔLᵈ > 0

                # TODO do this in a more sophisticated way
                n_sample = min(10, length(labormarket_struct.unemployed))

                # Make queue of job-seeking households
                Lˢ = sample(labormarket_struct.unemployed, n_sample, replace=false)

                # Only select households with a low enough reservation wage
                Lᵈ = sort(Lˢ, by = hh_id -> model[hh_id].wʳ)

                to_be_hired = Vector{Int}()
                ΔL = p.ΔLᵈ

                for hh_id in Lᵈ
                    
                    ΔL -= model[hh_id].L
                    push!(to_be_hired, hh_id)

                    if ΔL < 0
                        break
                    end
                end
                
                # Set offered wage to lowest requested wage that makes producer
                # meet the labor target
                p.wᴼ = model[to_be_hired[end]].wʳ

                # Hire selected workers
                for hh_id in to_be_hired

                    get_hired_hh!(model[hh_id], p.wᴼ, p_id)
                    hire_worker_p!(model[p_id], model[hh_id])                  

                end

                # Labor market aggregates are updated
                update_wage_level_p!(p, model)
                update_hiredworkers_lm!(labormarket_struct, to_be_hired)

                # If producer's labor demand is met, remove from seeking producers
                if p.ΔLᵈ < 0
                    filter!(p -> p ≠ p_id, hiring_producers)
                end

            end
        end
    end
end


function update_probs_lm!(
    labormarket_struct,
    model::ABM
    )

    n_unemployed = length(labormarket_struct.unemployed)
    n_longtermunemp = 0

    for hh_id in labormarket_struct.unemployed
        model[hh_id].T_unemp += 1
        if model[hh_id].T_unemp >= 2
            n_longtermunemp += 1
        end
    end

    labormarket_struct.P_UU = n_longtermunemp / n_unemployed
    labormarket_struct.P_HU = 1 - labormarket_struct.P_UU
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


function update_hiredworkers_lm!(
    labormarket_struct, 
    to_be_hired::Vector{Int}
    )

    # add newly employed workers to employed category
    append!(labormarket_struct.employed, to_be_hired)

    # remove employed workers from unemployed category
    filter!(hh -> hh ∉ to_be_hired, labormarket_struct.unemployed)
end


function update_firedworker_lm!(
    labormarket_struct,
    to_be_fired::Vector{Int}
    )

    # add newly unemployed workers to unemployed category
    append!(labormarket_struct.unemployed, to_be_fired)

    # remove unemployed workers from employed category
    filter!(hh -> hh ∉ to_be_fired, labormarket_struct.employed)

end