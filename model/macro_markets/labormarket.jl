@Base.kwdef mutable struct LaborMarket
    employed::Vector{Int} = []                # array of employed households
    unemployed::Vector{Int} = []              # array of unemployed households
    E::Float64 = 0.1                          # unemployment rate
    switch_rate::Float64 = 0.0                # rate by which employed employees switch employer
    n_rounds::Int = 3                         # number of rounds in matching process
    avg_T_unemp::Float64 = 0                  # average time periods of unemployment
end


"""
Gives all producers a share of the labor force after being initialized.
Households that are not assigned to a producer are declared unemployed.
"""
function spread_employees_lm!(
    labormarket_struct::LaborMarket,
    gov_struct::Government, 
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
        employees = all_hh[i:i+n_init_emp_cp-1]
        for hh_id in employees
            # w = max(gov_struct.w_min, model[hh_id].skill)
            w = 1.0
            set_employed_hh!(model[hh_id], w, cp_id)
            model[hh_id].w .= w
            hire_worker_p!(model[cp_id], model[hh_id])
            push!(labormarket_struct.employed, model[hh_id].id)
        end

        model[cp_id].employees = employees
        model[cp_id].L = sum(hh_id -> model[hh_id].L, employees)
        i += n_init_emp_cp

    end

    # Loop over kp, allocate households as employees
    for kp_id in all_kp
        employees = all_hh[i:i+n_init_emp_kp-1]
        for hh_id in employees
            w = max(gov_struct.w_min, model[hh_id].skill)
            set_employed_hh!(model[hh_id], w, kp_id)
            model[hh_id].w .= w
            hire_worker_p!(model[kp_id], model[hh_id])
            push!(labormarket_struct.employed, model[hh_id].id)
        end
        model[kp_id].employees = employees
        model[kp_id].L = sum(hh_id -> model[hh_id].L, employees)
        i += n_init_emp_kp
    end

    # Other households are pushed into unemployment
    for hh_id in all_hh[i:end]
        set_unemployed_hh!(model[hh_id])
        model[hh_id].w .= 0.0        
        push!(labormarket_struct.unemployed, model[hh_id].id)
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
    global_param::GlobalParam,
    gov_struct::Government,
    model::ABM,
    to
    )

    # get sets of firing and hiring producers
    firing_producers = Vector{Int}()
    hiring_producers = Vector{Int}()

    for p_id in all_p
        if model[p_id].ΔLᵈ > 0
            push!(hiring_producers, p_id)
        elseif model[p_id].ΔLᵈ < 0
            push!(firing_producers, p_id)
        end
    end

    # Let producers fire excess workers
    fire_workers_lm!(labormarket_struct, firing_producers, model)

    # Update wage parameters households
    for hh_id in all_hh
        update_sat_req_wage_hh!(model[hh_id], global_param.ϵ, gov_struct.UB, gov_struct.w_min)
    end

    if !global_param.fordist_lm
        # Find all employed households that want to change jobs
        employed_jobseekers = find_employed_jobseekers_lm(labormarket_struct.employed, global_param.ψ_E)
    else
        # In the Fordist regime, employed workers do not look for
        # a better paying job.
        employed_jobseekers = Vector{Int}()
    end

    # Labor market matching process
    @timeit to "matching" matching_lm(
        labormarket_struct, 
        employed_jobseekers, 
        hiring_producers,
        global_param.max_g_wᴼ, 
        model,
        to
    )

    # Update the unemployment rate
    update_unemploymentrate_lm!(labormarket_struct)
end


"""
Updates the unemployment rate in the economy
"""
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
    model::ABM,
    to
    )

    # Concatenate all unemployed and job-seeking employed households
    jobseeking_hh = vcat(employed_jobseekers, labormarket_struct.unemployed)
    hiring_producers_dict = Dict(p_id => model[p_id].ΔLᵈ for p_id in hiring_producers)

    # Track employed jobseekers that actually switch jobs
    n_jobswitchers = 0
    n_employed = length(labormarket_struct.employed)

    allowed_excess_L = 50

    for _ in 1:labormarket_struct.n_rounds

        # Loop over hiring producers producers
        for (p_id, ΔL) in hiring_producers_dict

            # Stop process if no unemployed left
            if length(jobseeking_hh) == 0
                return
            end

            # Make queue of job-seeking households
            n_sample = min(30, length(jobseeking_hh))
            @timeit to "sample" Lᵈ = sample(jobseeking_hh, n_sample, replace=false)
            @timeit to "sort" sort!(Lᵈ, by = hh_id -> model[hh_id].L * model[hh_id].skill, rev=true)

            to_be_hired = Vector{Int}()
            for hh_id in Lᵈ
                
                # Only hire workers if wage can be afforded
                # TODO: DESCRIBE IN MODEL
                if ((model[hh_id].wʳ / model[hh_id].skill <= model[p_id].wᴼ_max) &&
                    (model[hh_id].L * model[hh_id].skill <= ΔL + allowed_excess_L))
                    push!(to_be_hired, hh_id)
                    ΔL -= model[hh_id].L * model[hh_id].skill
                end

                # TODO: make this a constant
                # If producer's labor demand is met, remove from seeking producers
                if ΔL <= allowed_excess_L
                    filter!(p -> p ≠ p_id, hiring_producers)
                    break
                end
            end
            
            # Set offered wage to lowest requested wage that makes producer
            # meet the labor target

            unemployed_to_employed = Vector{Int}()

            if length(to_be_hired) > 0
                model[p_id].wᴼ = model[to_be_hired[end]].wʳ

                # Hire selected workers
                for hh_id in to_be_hired

                    # if model[hh_id].wʳ < 0.7
                    #     println(model[hh_id].wʳ)
                    # end

                    # If employed, change employer, otherwise get employed
                    if model[hh_id].employed
                        remove_worker_p!(model[model[hh_id].employer_id], model[hh_id])
                        # change_employer_hh!(model[hh_id], model[p_id].wᴼ, p_id)
                        change_employer_hh!(model[hh_id], model[hh_id].wʳ, p_id)
                        n_jobswitchers += 1
                    else
                        # set_employed_hh!(model[hh_id], model[p_id].wᴼ, p_id)
                        set_employed_hh!(model[hh_id], model[hh_id].wʳ, p_id)
                        push!(unemployed_to_employed, hh_id)
                    end
                    
                    hire_worker_p!(model[p_id], model[hh_id])
                    filter!(hh -> hh ≠ hh_id, jobseeking_hh)                  
                end

                # Labor market aggregates are updated
                update_w̄_p!(model[p_id], model)
                update_hiredworkers_lm!(labormarket_struct, unemployed_to_employed)

                hiring_producers_dict[p_id] = ΔL
            end
        end
    end

    # Updates the labor market's switching rate (use n employed from before matching)
    labormarket_struct.switch_rate = n_jobswitchers / n_employed
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
    filter!(hh_id -> hh_id ∉ to_be_fired, labormarket_struct.employed)

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