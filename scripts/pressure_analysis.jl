using DataFrames, FileIO, Plots

function read_Data(;p="lap", grad="step")
    go_back = "..\\..\\..\\..\\..\\"
    folder = "\\Software_Projects\\Swalbe.jl\\data\\Drop_coalescence_pressures\\"
    file = "$(p)_$(grad).jld2"
    data = load(go_back * folder * file) |> DataFrame
    return data
end

γ₀ = 1e-5

function step_gamma(;L=1024, γ=γ₀, perc=20, periodic=false, boundary_l=50)
    x = ones(L)
    for i in 1:L
        if i < L÷2
            x[i] = γ
        else
            x[i] = γ - (γ * perc / 100)
        end
    end
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

pressure_data = zeros(172, 3, 3, 2)
p_df = DataFrame()
plist = Float64[]
gnlist = String[]
pkind = String[]
tilist = Int64[]
for i in enumerate(["step", "tanh5", "tanh100"])
    gname = fill(i[2], 172)
    p_d = zeros(343)
    for n in enumerate(["lap", "disj"])
        pname = fill(n[2], 172)
        df_ = read_Data(p=n[2], grad=i[2])
    
        if i[2] == "step"
            gamma = step_gamma()
        else
            gamma = tanh_gamma(sl=parse(Int, i[2][5:end]))
        end
        for j in enumerate([10000,1000000,100000000])
            tname = fill(j[2], 172)
            p_d .= df_[!, Symbol("$(n[2])_$(j[2])")][512-171:512+171]
            calc .= p_d ./ (gamma[512-171:512+171] ./ 171)
            p_d .= calc
            pressure_data[:, i[1], j[1], n[1]] .= p_d[1:2:end]
            for l in 1:172
                push!(tilist, tname[l])
                push!(gnlist, gname[l])
                push!(plist, pressure_data[l, i[1], j[1], n[1]])
                push!(pkind, pname[l])
            end
        end
    end
end
p_df[!, "pressure"] = plist
p_df[!, "kind"] = pkind
p_df[!, "gamma"] = gnlist
p_df[!, "time"] = tilist

CSV.write("data\\pressure_data.csv", p_df)

# Figure for publication
function tau_ic(;ρ=1, r₀=171, γ=1e-5)
	return sqrt(ρ*r₀^3/γ)
end

