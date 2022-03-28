mutable struct Machine
    A :: Float64                # labor productivity machine
    freq :: Float64             # freq machine owned by cp
    age :: Float64              # age of machine
    p :: Float64                # price for which the machine was bought
end

"""
Initializes a machine struct.
"""
function initialize_machine(
    freq::Int,
    η::Int,
    p::Float64, 
    A::Float64
    )::Machine

    machine_struct = Machine(
        A,                      # A: labor productivity machine
        freq,                   # freq: freq machine owned by cp
        sample(0:η),            # age: age of machine
        p                       # p: price for which the machine was bought
    )
    return machine_struct
end


"""
Initializes a stock of machines for a new cp firm.
"""
function initialize_machine_stock(
    freq_per_machine::Int,
    n_machines_init::Int;
    η=0::Int,
    p=1.2::Float64,
    A=1.0::Float64,
    )::Vector{Machine}

    machines = Vector{Machine}()
    for _ in 1:n_machines_init
        machine_struct = initialize_machine(freq_per_machine, η, p, A)
        push!(machines, machine_struct)
    end

    return machines
end