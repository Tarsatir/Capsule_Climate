mutable struct LaborMarket
    employed :: Vector{Int}                     # array of employed households
    unemployed :: Vector{Int}                   # array of unemployed households
    E :: Float64                                # unemployment rate
    n_rounds :: Int                             # number of rounds in matching process
    avg_T_unemp :: Float64                      # average time periods of unemployment
end


function initialize_labormarket()
    labormarket_struct = LaborMarket(
        [],
        [],
        0.1,
        1,
        0
    )
    return labormarket_struct
end


"""
Gives all producers a share of the labor force after being initialized.
Households that are not assigned to a producer are declared unemployed.
"""
function spread_employees_lm!(
    labormarket_struct::LaborMarket, 
    all_hh::Vector{Int}, 
    all_cp::Vector{Int}, 
    all_kp::Vector{Int}, 
    n_init_emp_cp::Int, 
    n_init_emp_kp::Int,
    model::ABM
    )

    i = 1

    # Loop over cp, allocate households as employees
    for cp_id in all_cp
        cp = model[cp_id]
        employees = all_hh[i:i+n_init_emp_cp-1]
        for hh_id in employees
            hh = model[hh_id]
            set_employed_hh!(hh, 1.0, cp_id)
            hire_worker_p!(cp, hh)
            push!(labormarket_struct.employed, hh.id)
        end

        cp.employees = employees
        cp.L = sum(map(hh_id -> model[hh_id].L, employees))
        i += n_init_emp_cp

    end

    # Loop over kp, allocate households as employees
    for kp_id in all_kp
        kp = model[kp_id]
        employees = all_hh[i:i+n_init_emp_kp-1]
        for hh_id in employees
            hh = model[hh_id]
            set_employed_hh!(hh, 1.0, kp_id)
            hire_worker_p!(kp, hh)
            push!(labormarket_struct.employed, hh.id)
        end
        kp.employees = employees
        kp.L = sum(map(hh_id -> model[hh_id].L, employees))
        i += n_init_emp_kp
    end

    # Other households are pushed into unemployment
    for hh_id in all_hh[i:end]
        hh = model[hh_id]
        set_unemployed_hh!(hh)
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
    labormarket_struct::LaborMarket,
    all_hh::Vector{Int}, 
    all_p::Vector{Int}, 
    ϵ :: Float64,
    max_g_wᴼ :: Float64, 
    UB :: Float64,
    model::ABM
    )

    # get sets of firing and hiring producers
    firing_producers = Vector{Int}()
    hiring_producers = Vector{Int}()

    update_unemploymentrate_lm!(labormarket_struct)
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
    # println(sum(map(p_id -> length(model[p_id].employees), all_p)))
    fire_workers_lm!(labormarket_struct, firing_producers, model)
    # println("f2 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].employees), all_p)))

    # Update wage parameters households
    for hh_id in all_hh
        update_sat_req_wage_hh!(model[hh_id], ϵ, UB)
    end

    update_unemploymentrate_lm!(labormarket_struct)
    println("E 2: ", labormarket_struct.E)

    # Find all employed households that want to change jobs
    employed_jobseekers = find_employed_jobseekers_lm(labormarket_struct.employed, global_param.ψ_E)

    # labor market matching process
    # println("m1 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].employees), all_p)))
    matching_lm(
        labormarket_struct, 
        employed_jobseekers, 
        hiring_producers,
        max_g_wᴼ, 
        model
    )
    # println("m2 ", length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))
    # println(sum(map(p_id -> length(model[p_id].employees), all_p)))

    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

    # update the unemployment rate
    update_unemploymentrate_lm!(labormarket_struct)
    println("E 3: ", labormarket_struct.E)

    # println(length(labormarket_struct.employed), " ", length(labormarket_struct.unemployed))

end


function update_unemploymentrate_lm!(
    labormarket_struct::LaborMarket
    )
    labormarket_struct.E = (length(labormarket_struct.unemployed) / (length(labormarket_struct.employed) + length(labormarket_struct.unemployed)))
end


"""
Employers select workers to fire, labor market aggregates are changed and employees 
are fired.
"""
function fire_workers_lm!(
    labormarket_struct::LaborMarket, 
    firing_producers::Vector{Int}, 
    model::ABM
    )

    all_fired_workers = Vector{Int}()

    # Choose who gets fired
    for p_id in firing_producers
        fired_workers = fire_excess_workers_p!(model[p_id], model)
        # println(fired_workers)
        append!(all_fired_workers, fired_workers)
    end

    # Update employed and unemployed lists
    update_firedworker_lm!(labormarket_struct, all_fired_workers)

    # Change employment status for households
    for hh_id in all_fired_workers
        set_unemployed_hh!(model[hh_id])
    end
end


"""
Matches job-seeking households with employee-seeking firms
"""
function matching_lm(
    labormarket_struct::LaborMarket,
    employed_jobseekers::Vector{Int}, 
    hiring_producers::Vector{Int},
    max_g_wᴼ::Float64,
    model::ABM
    )

    # Concatenate all unemployed and job-seeking employed households
    jobseeking_hh = vcat(employed_jobseekers, labormarket_struct.unemployed)
    hiring_producers_dict = Dict(p_id => model[p_id].ΔLᵈ for p_id in hiring_producers)

    for i_round in 1:labormarket_struct.n_rounds

        # Loop over hiring producers producers
        for (p_id, ΔL) in hiring_producers_dict

            # Stop process if no unemployed left
            if length(jobseeking_hh) == 0
                return
            end

            # TODO do this in a more sophisticated way
            n_sample = min(10, length(jobseeking_hh))

            # Make queue of job-seeking households
            Lˢ = sample(jobseeking_hh, n_sample, replace=false)

            # Only select households with a low enough reservation wage
            Lᵈ = sort(Lˢ, by = hh_id -> model[hh_id].wʳ)

            to_be_hired = Vector{Int}()
            for hh_id in Lᵈ
                
                # Only hire workers if wage can be afforded
                if model[hh_id].wʳ <= model[p_id].wᴼ * (1 + max_g_wᴼ)
                    push!(to_be_hired, hh_id)
                    ΔL -= model[hh_id].L
                end

                # TODO: make this a constant
                # If producer's labor demand is met, remove from seeking producers
                if ΔL <= 10
                    filter!(p -> p ≠ p_id, hiring_producers)
                    break
                end
            end
            
            # Set offered wage to lowest requested wage that makes producer
            # meet the labor target

            if length(to_be_hired) > 0
                model[p_id].wᴼ = model[to_be_hired[end]].wʳ

                # Hire selected workers
                for hh_id in to_be_hired

                    # If employed, change employer, otherwise get employed
                    if model[hh_id].employed
                        remove_worker_p!(model[model[hh_id].employer_id], model[hh_id])
                        change_employer_hh!(model[hh_id], model[p_id].wᴼ, p_id)
                    else
                        set_employed_hh!(model[hh_id], model[p_id].wᴼ, p_id)
                    end
                    
                    hire_worker_p!(model[p_id], model[hh_id])
                    filter!(hh -> hh ≠ hh_id, jobseeking_hh)                  
                end

                # Labor market aggregates are updated
                update_w̄_p!(model[p_id], model)
                update_hiredworkers_lm!(labormarket_struct, to_be_hired)

                hiring_producers_dict[p_id] = ΔL
            end
        end
    end
end


# function update_avg_T_unemp_lm(
#     labormarket_struct, 
#     model::ABM
#     )

#     if length(labormarket_struct.unemployed) > 0
#         avg_T_unemp = mean(map(hh_id -> model[hh_id].T_unemp, labormarket_struct.unemployed))
#         labormarket_struct.avg_T_unemp = avg_T_unemp
#     else
#         labormarket_struct.avg_T_unemp = 0
#     end
# end


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
    labormarket_struct::LaborMarket,
    to_be_fired::Vector{Int}
    )

    # add newly unemployed workers to unemployed category
    append!(labormarket_struct.unemployed, to_be_fired)

    # remove unemployed workers from employed category
    filter!(hh -> hh ∉ to_be_fired, labormarket_struct.employed)

end


"""
Finds employed workers that want to look for another job.
"""
function find_employed_jobseekers_lm(
    employed::Vector{Int},
    ψ_E::Float64
    )::Vector{Int}

    n = floor(Int, length(employed) * ψ_E)
    employed_jobseekers = sample(employed, n)
    return employed_jobseekers
end