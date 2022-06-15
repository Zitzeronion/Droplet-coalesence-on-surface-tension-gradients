using DataFrames, FileIO, Plots

function read_Data(;p="lap", grad="step")
    go_back = "..\\..\\..\\..\\..\\"
    folder = "\\Software_Projects\\Swalbe.jl\\data\\Drop_coalescence_pressures\\"
    file = "$(p)_$(grad).jld2"
    data = load(go_back * folder * file) |> DataFrame
    return data
end

df = read_Data()
t = 1000000
plot(df[!, Symbol("lap_$(t)"], label="Δh @ $(t)")