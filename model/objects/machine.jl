mutable struct Machine
    A :: Float64                # labor productivity machine
    freq :: Float64             # freq machine owned by cp
    age :: Float64              # age of machine
    p :: Float64                # price for which the machine was bought
end

"""
Initializes machine struct
"""
function initialize_machine(
    freq::Float64;
    η=0::Int,
    p=1.0::Float64, 
    A=1.0::Float64
    )

    machine_struct = Machine(
        A,                      # A: labor productivity machine
        freq,                   # freq: freq machine owned by cp
        sample(0:η),            # age: age of machine
        p                       # p: price for which the machine was bought
    )
    return machine_struct
end


function initialize_machine_stock(
    tot_freq_machines::Float64,
    p::Float64,
    A::Float64,
    n_machines_init::Int
    )

    machines = Vector{Machine}()
    for i in 1:n_machines_init
        freq = tot_freq_machines/n_machines_init
        machine_struct = initialize_machine(freq; η=0, p=p, A=A)
        push!(machines, machine_struct)
    end

    return machines
end