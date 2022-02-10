mutable struct Balance
    # assets
    N :: Float64            # inventories
    K :: Float64            # capital
    NW :: Float64           # liquid assets

    # liabilities
    Deb :: Float64          # debt
    EQ :: Float64           # equity
end

mutable struct GovCurrentAccount
    # Revenues
    Rev_τᴵ :: Vector{Float64}            # hist revenues of income tax
    Rev_τˢ :: Vector{Float64}            # hist revenues of sales tax
    Rev_τᴾ :: Vector{Float64}            # hist revenues of profit tax
    Rev_τᴱ :: Vector{Float64}            # hist revenues of energy tax
    Rev_τᶜ :: Vector{Float64}            # hist revenues of emission tax

    # Expenditures
    Exp_UB :: Vector{Float64}            # hist expenditures on unemployment benefits
    Exp_Sub :: Vector{Float64}           # hist expenditures on subsidies

end

# function initialize_govcurrentaccount()
#     gov_curr_acc = GovCurrentAccount(
#         [],
#         [],
#         [],
#         [],
#         [],
#         [],
#         []
#     )
#     return gov_curr_acc
# end

mutable struct Machine
    A :: Float64                # labor productivity machine
    c :: Float64                # cost to produce machine
    freq :: Float64             # freq machine owned by cp
    age :: Float64              # age of machine
end

function initialize_machine()
    machine_struct = Machine(
        1,                      # A: labor productivity machine
        0,                      # c: cost to produce machine
        40,                     # freq: freq machine owned by cp
        0                       # age: age of machine
    )
    return machine_struct
end