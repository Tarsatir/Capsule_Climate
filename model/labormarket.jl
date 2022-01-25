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
        0.1,
        1
    )
    return labormarket_struct
end


"""
Defines the whole labor market process

- Lets producers fire excess workers
- Matches unemployed workers with producers looking for workers
- Updates unemployment rate

"""
function labormarket_process!(labormarket_struct, all_cp, all_kp)

    # get sets of firing and hiring producers
    firing_producers = []
    hiring_producers = []

    # update_unemploymentrate_lm(labormarket_struct)
    # println(labormarket_struct.E)

    for p in vcat(all_cp, all_kp)
        if p.ΔLᵈ > 0
            push!(hiring_producers, p)
        elseif p.ΔLᵈ < -100
            push!(firing_producers, p)
        end
    end

    # let producers fire excess workers
    fire_workers!(labormarket_struct, firing_producers)

    # update_unemploymentrate_lm(labormarket_struct)
    # println(labormarket_struct.E)

    # labor market matching process
    matching_lm(labormarket_struct, all_cp, all_kp)

    # update the unemployment rate
    update_unemploymentrate_lm(labormarket_struct)
    # println(labormarket_struct.E)

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
        append!(all_fired_workers, fire_excess_workers!(p))
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

    # get producers that want to hire workers
    all_p = []
    for p in vcat(all_cp, all_kp)
        if p.ΔLᵈ > 0
            push!(all_p, p)
        end
    end

    for n_round in labormarket_struct.n_rounds

        # loop over producers
        shuffle!(all_p)
        for p in all_p

            # TODO do this in a more sophisticated way
            n_sample = min(10, length(Lᵃ))

            # make queue of job-seeking households
            Lˢ = sample(Lᵃ, n_sample, replace=false)

            # only select households with a low enough reservation wage
            Lᵈ = []
            for l in Lˢ
                if l.wʳ <= p.wᴼ
                    push!(Lᵈ, l)
                end
            end

            # hire workers until demand is met or no more workers available
            while p.ΔLᵈ > 0 && length(Lᵈ) > 0

                # hire worker
                l = Lᵈ[1]
                get_hired_hh!(l, p)
                hire_worker_p!(p, l)

                # delete household from seeking lists
                filter!(w -> w != l, Lᵈ)
                filter!(w -> w != l, Lᵃ)

            end
        end
    end
end