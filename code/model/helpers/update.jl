"""
    shift_and_append!(ts::Vector{Union{Float64, Int64}}, neval::Union{Float64, Int64})

Shifts all elements in the passed array to the left and adds the new
    value to the end.
"""
function shift_and_append!(
    ts::Union{Vector{Float64}, Vector{Int64}},
    newval::Union{Float64, Int64}
    )

    # Shift values
    ts[1:end-1] = ts[2:end]

    # Add new value
    ts[end] = newval
end