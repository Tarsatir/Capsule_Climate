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
    freq::Float64,
    η::Int=0,
    p::Float64=1.0, 
    A::Float64=1.0
    )

    machine_struct = Machine(
        A,                      # A: labor productivity machine
        freq,                   # freq: freq machine owned by cp
        sample(0:η),            # age: age of machine
        p                       # p: price for which the machine was bought
    )
    return machine_struct
end