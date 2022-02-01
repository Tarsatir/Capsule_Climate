"""
File used to write simulation results to data files
"""

using DataFrames
using CSV

"""
Saves macro variables of interest to csv

Receives:
    macro_struct: mut struct with macro variables of interest
"""
function save_macro_data(macro_struct)

    df = DataFrame(GDP=macro_struct.GDP,
                   UR=macro_struct.U)
    CSV.write("results/result_data/first.csv", df)
end