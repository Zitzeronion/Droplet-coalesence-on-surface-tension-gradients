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

# Figure for publication
function tau_ic(;ρ=1, r₀=171, γ=1e-5)
	return sqrt(ρ*r₀^3/γ)
end

p_list = Float64[]
gamma_l = String[]
kind_l = String[]
time_l = Int64[]
time_n = Float64[]

for i in enumerate(["step", "tanh5", "tanh100"])
    df_l = read_Data(p="lap", grad=i[2])
    df_pi = read_Data(p="disj", grad=i[2])
    if i[2] == "step"
        gamma = step_gamma()
    else
        gamma = tanh_gamma(sl=parse(Int, i[2][5:end]))
    end
    
    # c_ = [palette(:default)[1], palette(:default)[2], palette(:default)[3]]
    for j in enumerate([100000, 1000000, 10000000])
        for k in [df_l, df_pi]
            if k == df_l
                kind = "lap"
            else
                kind = "disj"
            end
            data[:, 1] .= k[!, Symbol("$(kind)_$(j[2])")][512-171:512+171]
            # data[:, 2] .= df_pi[!, Symbol("disj_$(j[2])")][512-171:512+171]
            calc .= data[:, 1] ./ (gamma[512-171:512+171] ./ 171)
            data[:, 1] .= calc
            # calc .= data[:, 2] ./ (gamma[512-171:512+171] ./ 171)
            # data[:, 2] .= calc
            for u in data[:, 1]
                push!(p_list, u)
                push!(gamma_l, i[2])
                push!(kind_l, kind)
                push!(time_l, j[2])
                push!(time_n, j[2]/tau_ic())
            end
        end
        #=
        xax = collect(-171:171) ./ 171
        plot!(xax, data[1:2:end, 1], 
            label="∂h/∂ₓ",
            l=(3, :solid, c_[j[1]]), 
            xlabel="x/R_0",
            ylabel="p/P") 
        plot!(xax, data[1:2:end, 2], l=(3, :dash, c_[j[1]]), label="Π")  
        =#
    end
    # savefig(p_dist, "figures\\pressures_dists_$(i[2]).svg")
end

df = DataFrame()
df[!, "P_value"] = p_list
df[!, "gamma"] = gamma_l
df[!, "kind"] = kind_l
df[!, "time"] = time_l
df[!, "time_n"] = time_n
# savefig(p_dist, "figures\\pressures_dists.svg")

function pressures(df; gamma="step", kind="lap", time=100000)
    return df[(df.gamma .== gamma) .& (df.kind .== kind) .& (df.time .== time), :P_value] 
end
R = 171
plot(collect(-R:R) ./ R, pressures(df, kind="disj", time=10000000))