@Base.kwdef mutable struct LaborMarket
    employed_hh::Vector{Int} = Int64[]        # array of employed households
    unemployed::Vector{Int} = Int64[]         # array of unemployed households
    jobseeking_hh::Vector{Int} = Int64[]      # array of jobseeking households

    hiring_producers::Vector{Int} = Int64[]   # array of hiring producers
    firing_producers::Vector{Int} = Int64[]   # array of firing producers 

    L_demanded::Float64 = 0.0                 # total labor demanded
    L_hired::Float64 = 0.0                    # total labor hired
    L_fired::Float64 = 0.0                    # total labor fired
    E::Float64 = 0.0                          # unemployment rate
    switch_rate::Float64 = 0.0                # rate by which employed employees switch employer
end


"""
Gives all producers a share of the labor force after being initialized.
Households that are not assigned to a producer are declared unemployed.
"""
function spread_employees_lm!(
    labormarket::LaborMarket,
    government::Government, 
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
            w = max(government.w_min, model[hh_id].skill)
            set_employed_hh!(model[hh_id], w, cp_id)
            model[hh_id].w .= w
            hire_worker_p!(model[cp_id], model[hh_id])
            push!(labormarket.employed_hh, hh_id)
        end

        model[cp_id].employees = employees
        model[cp_id].L = length(employees) > 0 ? sum(hh_id -> model[hh_id].L * model[hh_id].skill, employees) : 0.0
        i += n_init_emp_cp

    end

    # Loop over kp, allocate households as employees
    for kp_id in all_kp
        employees = all_hh[i:i+n_init_emp_kp-1]

        for hh_id in employees
            w = max(government.w_min, model[hh_id].skill)
            set_employed_hh!(model[hh_id], w, kp_id)
            model[hh_id].w .= w
            hire_worker_p!(model[kp_id], model[hh_id])
            push!(labormarket.employed_hh, hh_id)
        end
        model[kp_id].employees = employees
        model[kp_id].L = length(employees) > 0 ? sum(hh_id -> model[hh_id].L * model[hh_id].skill, employees) : 0.0
        i += n_init_emp_kp
    end

    # Other households are pushed into unemployment
    for hh_id in all_hh[i:end]
        set_unemployed_hh!(model[hh_id])
        model[hh_id].w .= 0.0        
        push!(labormarket.unemployed, model[hh_id].id)
    end

    update_unemploymentrate_lm!(labormarket)
end


"""
Defines the whole labor market process

- Lets producers fire excess workers
- Matches unemployed workers with producers looking for workers
- Updates unemployment rate

"""
function labormarket_process!(
    labormarket::LaborMarket,
    all_hh::Vector{Int}, 
    all_p::Vector{Int}, 
    globalparam::GlobalParam,
    government::Government,
    t::Int,
    model::ABM,
    to
    )

    # Update which producers want to hire or fire workers
    update_hiring_firing_producers(labormarket, all_p, model)

    # Let producers fire excess workers
    fire_workers_lm!(labormarket, model)

    # Update wage parameters households
    for hh_id in all_hh
        # update_sat_req_wage_hh!(
        #     model[hh_id], 
        #     globalparam.ϵ, 
        #     globalparam.ω, 
        #     government.w_min
        # )

        update_sat_req_wage_hh!(
            model[hh_id], 
            globalparam.ϵ_w,
            government.w_min
        )
    end

    # Find all employed households that want to change jobs
    employed_jobseekers = find_employed_jobseekers_lm(labormarket.employed_hh, globalparam.ψ_E)

    # Update jobseeking households
    labormarket.jobseeking_hh = vcat(employed_jobseekers, labormarket.unemployed)

    # Labor market matching process
    @timeit to "matching" matching_lm(
        labormarket,
        all_p, 
        model,
    )

    # Update the unemployment rate
    update_unemploymentrate_lm!(labormarket)

    # Round labor amounts to avoid floating point errors
    for p_id in all_p
        model[p_id].L = length(model[p_id].employees) == 0 ? 0.0 : model[p_id].L
        update_L!(model[p_id], model)
    end
end


"""
    update_hiring_firing_producers(labormarket::LaborMarket, all_p::Vector{Int}, model::ABM)

Updates the arrays of hiring and firing producers in the labor market struct.    
"""
function update_hiring_firing_producers(
    labormarket::LaborMarket,
    all_p::Vector{Int},
    model::ABM
    )

    # Get sets of firing and hiring producers
    labormarket.firing_producers = Int64[]
    labormarket.hiring_producers = Int64[]

    # TODO: write this in a nicer form
    for p_id in all_p
        if model[p_id].ΔLᵈ > 0
            push!(labormarket.hiring_producers, p_id)
        elseif model[p_id].ΔLᵈ < 0
            push!(labormarket.firing_producers, p_id)
        end
    end

end


"""
Updates the unemployment rate in the economy
"""
function update_unemploymentrate_lm!(
    labormarket::LaborMarket
    )

    labormarket.E = (length(labormarket.unemployed) / (length(labormarket.employed_hh) + length(labormarket.unemployed)))
end


"""
Employers select workers to fire, labor market aggregates are changed and employees 
are fired.
"""
function fire_workers_lm!(
    labormarket::LaborMarket, 
    model::ABM
    )

    all_fired_workers = Int64[]

    # Choose who gets fired
    for p_id in labormarket.firing_producers
        fired_workers = fire_excess_workers_p!(model[p_id], model)
        append!(all_fired_workers, fired_workers)
    end

    # Update employed and unemployed lists
    update_firedworker_lm!(labormarket, all_fired_workers)

    # Change employment status for households
    for hh_id in all_fired_workers
        set_unemployed_hh!(model[hh_id])
    end

    labormarket.L_fired = length(all_fired_workers) > 0 ? sum(hh_id -> model[hh_id].L * model[hh_id].skill, all_fired_workers) : 0.0
end


"""
Matches job-seeking households with employee-seeking firms
"""
function matching_lm(
    labormarket::LaborMarket,
    all_p::Vector{Int},
    model::ABM,
    )
    
    # jobseeking_weights = map(hh_id -> model[hh_id].skill, labormarket.jobseeking_hh)

    # Track employed jobseekers that actually switch jobs
    n_jobswitchers = 0
    n_employed = length(labormarket.employed_hh)

    allowed_excess_L = 100

    Lᵈ = Vector{Int}(undef, 100)

    labormarket.L_demanded = length(labormarket.hiring_producers) > 0 ? sum(p_id -> model[p_id].ΔLᵈ, labormarket.hiring_producers) : 0.0
    labormarket.L_hired = 0.0

    # update_hiring_firing_producers(labormarket, all_p, model)
    sort!(labormarket.hiring_producers, by=p_id->model[p_id].ΔLᵈ, rev=true)
    hiring_producers_dict = Dict(p_id => model[p_id].ΔLᵈ for p_id in labormarket.hiring_producers)

    # Loop over hiring producers producers
    for (p_id, ΔL) in hiring_producers_dict
        # println(ΔL)

        demanded_labor = ΔL

        # Stop process if no unemployed left
        if length(labormarket.jobseeking_hh) == 0
            return
        elseif length(labormarket.jobseeking_hh) < length(Lᵈ)
            Lᵈ = Vector{Int}(undef, length(labormarket.jobseeking_hh))
        end

        # Make queue of job-seeking households
        # n_sample = min(30, length(jobseeking_hh))
        # @timeit to "sample" Lᵈ .= sample(jobseeking_hh, n_sample, replace=false)
        StatsBase.sample!(labormarket.jobseeking_hh, Lᵈ; replace=false)
        # sort!(Lᵈ, by = hh_id -> model[hh_id].skill, rev=true)
        # w = map(hh_id -> model[hh_id].skill, labormarket.jobseeking_hh)
        # StatsBase.sample!(Lᵈ, Lᵈ, weights(w); replace=false)

        to_be_hired = Vector{Int}()

        for hh_id in Lᵈ
            # Only hire workers if wage can be afforded
            # TODO: DESCRIBE IN MODEL
            if ((model[hh_id].wʳ / model[hh_id].skill <= model[p_id].wᴼ_max) &&
                (model[hh_id].L * model[hh_id].skill <= demanded_labor + allowed_excess_L))
                push!(to_be_hired, hh_id)
                demanded_labor -= model[hh_id].L * model[hh_id].skill
                labormarket.L_hired += model[hh_id].L * model[hh_id].skill
            end

            # TODO: make this a constant
            # If producer's labor demand is met, remove from seeking producers
            if demanded_labor <= allowed_excess_L
                filter!(p -> p ≠ p_id, labormarket.hiring_producers)
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

                # If employed, change employer, otherwise get employed
                if model[hh_id].employed
                    remove_worker_p!(model[model[hh_id].employer_id], model[hh_id])
                    change_employer_hh!(model[hh_id], model[hh_id].wʳ, p_id)
                    n_jobswitchers += 1
                else
                    set_employed_hh!(model[hh_id], model[hh_id].wʳ, p_id)
                    push!(unemployed_to_employed, hh_id)
                end
                
                hire_worker_p!(model[p_id], model[hh_id])
                filter!(hh -> hh ≠ hh_id, labormarket.jobseeking_hh)                  
            end

            # Labor market aggregates are updated
            update_w̄_p!(model[p_id], model)
            update_hiredworkers_lm!(labormarket, unemployed_to_employed)

            hiring_producers_dict[p_id] = demanded_labor
        end
    end

    # Updates the labor market's switching rate (use n employed from before matching)
    labormarket.switch_rate = n_jobswitchers / n_employed
end


function update_hiredworkers_lm!(
    labormarket::LaborMarket, 
    to_be_hired::Vector{Int}
    )

    # Add newly employed workers to employed category
    append!(labormarket.employed_hh, to_be_hired)

    # Remove employed workers from unemployed category
    filter!(hh_id -> hh_id ∉ to_be_hired, labormarket.unemployed)
end


function update_firedworker_lm!(
    labormarket::LaborMarket,
    to_be_fired::Vector{Int}
    )

    # Add newly unemployed workers to unemployed category
    append!(labormarket.unemployed, to_be_fired)

    # Remove unemployed workers from employed category
    filter!(hh_id -> hh_id ∉ to_be_fired, labormarket.employed_hh)
end


"""
Finds employed workers that want to look for another job.
"""
function find_employed_jobseekers_lm(
    employed_hh::Vector{Int},
    ψ_E::Float64
    )::Vector{Int}

    n = floor(Int64, length(employed_hh) * ψ_E)
    employed_jobseekers = sample(employed_hh, n)
    return employed_jobseekers
end