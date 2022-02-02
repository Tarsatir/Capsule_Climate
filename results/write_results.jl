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
                   GDP_I=macro_struct.GDP_I,
                   GDP_cp=macro_struct.GDP_Π_cp,
                   GDP_kp=macro_struct.GDP_Π_kp,
                   UR=macro_struct.U,
                   s_avg=macro_struct.s̄_avg,
                   s_std=macro_struct.s̄_std,
                   w_avg=macro_struct.w̄_avg,
                   w_std=macro_struct.w̄_std)
    CSV.write("results/result_data/first.csv", df)
end