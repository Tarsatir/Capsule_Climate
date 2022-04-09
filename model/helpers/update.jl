"""
Shifts all elements in the passed array to the left and adds the new
    value to the end.
"""
function shift_and_append!(
    ts::Vector{Float64},
    newval::Float64
    )

    # Shift values
    ts[1:end-1] = ts[2:end]

    # Add new value
    ts[end] = newval
end