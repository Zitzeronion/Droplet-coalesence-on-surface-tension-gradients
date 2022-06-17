using DataFrames, FileIO, Plots

function read_Data(;p="lap", grad="step")
    go_back = "..\\..\\..\\..\\..\\"
    folder = "\\Software_Projects\\Swalbe.jl\\data\\Drop_coalescence_pressures\\"
    file = "$(p)_$(grad).jld2"
    data = load(go_back * folder * file) |> DataFrame
    return data
end

γ₀ = 1e-5

function tanh_gamma(;L=1024, γ=γ₀, perc=20, sl=1, periodic=false, boundary_l=50)
    l = collect(1.0:L)
    function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	x = ones(L)
	x .= γ .* smooth(l, L, sl) .+ (1 .- smooth(l, L, sl)) .* (γ - (γ * perc / 100)) 
	# Linear interpolation between surface tensions at the boundary
	if periodic
		for i in enumerate(boundary_l:-1:1)
			x[i[1]] = γ - i[2]/boundary_l * Δ/2
		end
		for i in enumerate(L-boundary_l:L)
			x[i[2]] = γ - Δ + i[1]/boundary_l * Δ/2
		end
	end
	return x
end

data = zeros(343, 2)
calc = zeros(343, 1)

for i in ["tanh5", "tanh100"]
    df_l = read_Data(p="lap", grad=i)
    df_pi = read_Data(p="disj", grad=i)
    gamma = tanh_gamma(sl=parse(Int, i[5:end]))
    drops = @animate for i in 10000:100000:100000000
        data[:, 1] .= df_l[!, Symbol("lap_$(i)")][512-171:512+171]
        data[:, 2] .= df_pi[!, Symbol("disj_$(i)")][512-171:512+171]
        calc .= data[:, 1] ./ (gamma[512-171:512+171] ./ 171)
        data[:, 1] .= calc
        calc .= data[:, 2] ./ (gamma[512-171:512+171] ./ 171)
        data[:, 2] .= calc 
        plot(data[1:2:end, :], label=["lap" "Pi"])
    end
    gif(drops, "figures\\pressures_$(i).gif")
end
