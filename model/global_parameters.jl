@with_kw mutable struct GlobalParam
    # Determine technical innovation process
    ν::Float64 = 0.05               # R&D inv propensity (note this a greek 'nu')
    νₑ::Float64 = 0.01              # R&D inv propensity energy producer (note this a greek 'nu')
    ξ::Float64 = 0.3                # R&D allocation to IN
    ξₑ::Float64 = 0.4               # R&D allocation to green tech for energy producer
    ζ::Float64 = 0.3                # firm search capabilities
    ζ_ge::Float64 = 0.3             # ep search capabilities for green tech
    ζ_de::Float64 = 0.3             # ep search capabilities for dirty tech
    α1::Float64 = 3.0               # 1st beta dist param for IN
    β1::Float64 = 3.0               # 2nd beta dist param for IN

    κ_upper::Float64 = 0.005        # 2nd beta dist support
    κ_lower::Float64 = -κ_upper     # 1st beta dist support
    
    γ::Float64 = 0.5                # new custommer sample parameter
    μ1::Float64 = 0.2               # kp markup rule
    v::Float64 = 0.04               # 
    r::Float64 = 0.0                # interest rate
    ι::Float64 = 0.2                # desired inventories
    b::Int64 = 9                    # payback period
    bₑ::Int64 = 30                  # payback period energy producer
    η::Int64 = 60                   # physical scrapping age
    ηₑ::Int64 = 240                 # physical scrapping age energy producer
    Λ::Float64 = 2.0                # max debt/sales ratio
    update_period::Int64 = 3        # time period after which cp update prod plans

    # Determine entrant composition
    φ1::Float64 = 0.1               # 1st Uniform dist support, cp entrant cap
    φ2::Float64 = 0.9               # 2nd Uniform dist support, cp entrant cap
    φ3::Float64 = 0.1               # 1st Uniform dist support, cp entrant liq
    φ4::Float64 = 0.9               # 2nd Uniform dist support, cp entrant liq
    α2::Float64 = 2.0               # 1st beta dist param for kp entrant
    β2::Float64 = 4.0               # 2nd beta dist param for kp entrant
    φ5::Float64 = -0.04             # 1st Beta dist support for kp entrant tech
    φ6::Float64 = 0.02              # 2nd Beta dist support for kp entrant tech

    cu::Float64 = 0.75              # capacity utilization for cp
    max_NW_ratio::Float64 = 0.5     # maximum ratio p can have monthly expenses in NW
    ϵ_w::Float64 = 0.03             # minimum desired wage increase rate
    ϵ_μ::Float64 = 0.01             # upper limit of markup shock
    Kg_max::Float64 = 0.5           # maximum capital growth rate

    # Determine expectation updating cp
    ω::Float64 = 0.8                # memory parameter adaptive updating rules
    λ::Float64 = 0.7                # parameter for labor demand smoothing
    p_rigid_time::Int64 = 1         # number of periods producer's price remains stable

    # Determine household consumption
    α_cp::Float64 = 0.75            # parameter controlling APC of consumers

    # Deterime extend of progresivity of government spending
    prog::Float64 = -0.5

    # Determine household switching
    ψ_E::Float64 = 0.15             # chance of employed worker looking for a better paying job
    ψ_Q::Float64 = 0.05             # chance of household switching away from cp when demand constrained
    ψ_P::Float64 = 0.05             # chance of household switching to cp with better price

    freq_per_machine::Int64 = 25    # capital units per machine
    freq_per_powerplant::Int64 = 10_000 # capital units per instance

    p_f::Float64 = 0.2              # price of fossil fuels

    n_cons_market_days::Int64 = 4   # number of days in the consumer market process

    t_warmup::Int64 = 300           # time period warmup of the model
    t_wait::Int64 = 4               # number of time periods producers are not allowed to go bankrupt

    changedparams_ofat::Union{Nothing, Dict} # Parameters that are changed at the end of the warmup period
end


function initialize_global_params(
    changedparams::Union{Nothing, Dict},
    changedparams_ofat::Union{Nothing, Dict}
    )

    globalparam = GlobalParam(changedparams_ofat=changedparams_ofat)

    # Change parameters if needed before returning.
    if changedparams ≠ nothing
        for (key, new_param) in changedparams
            setproperty!(globalparam, Symbol(key), new_param)
        end

        if haskey(changedparams, "κ_upper")
            setproperty!(globalparam, Symbol("κ_lower"), -changedparams["κ_upper"])
        end
    end

    return globalparam
end


"""
Checks if parameters are changed during the simulation. If they are, checks if time of
    introduction has been reached. If time reached, changes parameter value.
"""
function check_changed_ofatparams(
    globalparam::GlobalParam,
    t::Int64
    )

    if globalparam.changedparams_ofat ≠ nothing
        for (paramtype, (newparamval, oldparamval, t_introduction, t_duration)) in globalparam.changedparams_ofat
            if t == t_introduction
                # Set old value to current value of parameter
                globalparam.changedparams_ofat[paramtype][2] = getfield(globalparam, paramtype)

                # Set parameter to new value
                setproperty!(globalparam, paramtype, newparamval)
            elseif t == t_introduction + t_duration
                # Set parameter value back to old value
                setproperty!(globalparam, paramtype, oldparamval)
            end
        end
    end
end