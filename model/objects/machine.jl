@with_kw mutable struct Machine
    A_LP::Float64             # Labor productivity machine
    A_EE::Float64             # Energy efficiency machine
    A_EF::Float64             # Environmental friendliness machine
    freq::Float64             # Freq machine owned by cp
    age::Float64              # Age of machine
    p::Float64                # Price at which machine was bought
    cop::Float64              # Current cost of production of machine
end


"""
Initializes a machine struct.
"""
function initialize_machine(
    freq::Int,
    η::Int,
    p::Float64, 
    A_LP::Float64,
    A_EE::Float64,
    A_EF::Float64
    )::Machine

    machine_struct = Machine(
        A_LP=A_LP,
        A_EE=A_EE,
        A_EF=A_EF,
        freq=freq,
        age=sample(0:η),
        p=p,
        cop=0.0        
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
    A_LP=1.0::Float64,
    A_EE=1.0::Float64,
    A_EF=1.0::Float64
    )::Vector{Machine}

    # machines = Vector{Machine}()
    # for _ in 1:n_machines_init
    #     machine_struct = initialize_machine(freq_per_machine, η, p, A_LP, A_EE, A_EF)
    #     push!(machines, machine_struct)
    # end

    machines = [initialize_machine(freq_per_machine, η, p, A_LP, A_EE, A_EF) for _ in 1:n_machines_init]

    return machines
end


function update_cop_machine!(
    machine::Machine,
    cp::AbstractAgent,
    government,
    ep::AbstractAgent,
    t::Int64
)

    machine.cop = cop(cp.w̄[end], machine.A_LP, government.τᴱ, ep.pₑ[t], machine.A_EE, government.τᶜ, machine.A_EF)
end