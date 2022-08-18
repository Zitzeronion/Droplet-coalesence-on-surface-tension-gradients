### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 50901b90-eb08-11ec-3087-81a949a22d5d
using CSV, DataFrames, DataFramesMeta, Plots, PlutoUI, HTTP, FileIO, JLD2

# ╔═╡ 7ed11bdd-2f22-44a5-94c2-8ef0d841769e
md"# Flows and droplets

Here we want to take a look at the evolution of the two maxima.
If a droplet gets smaller than the other one should increase in size.
Measuring this is not as trivial as trivial as one expects.
Because the geometry is not a circular section all the time.
To make things simple we are going to use the metric

```math
	\Delta h = h_d^l - h_d^r
```

where $h_d$ is the maximum height of the droplets.
At $t=0$ we set $h_d^l = h_d^r$, but the right side $h_d^r$ is place on the region with lower surface tension.
We expect that $h_d^r$ will decrease with time. 
"

# ╔═╡ cb1d2d4f-5f05-4fe9-93db-973c9113cf67
md"First we read the data, simply clone the repo and load the .csv file."

# ╔═╡ ef661fd9-e3c2-4fd3-a36b-75a43222155f
#coalescene_data = CSV.read("..\\data\\coalescence_data_slip_12_gamma0_1e-5_sim.csv", DataFrame)
coalescene_data = CSV.read(HTTP.get("https://raw.githubusercontent.com/Zitzeronion/TimeDepWettabilityPaper/Drop_coalescence/Data_CSV/coalescence_data_slip_12_gamma0_1e-5_sim.csv").body, DataFrame)

# ╔═╡ 0671ef84-5f84-4967-bf29-b3f9e83f9f08
#coalescene_data2 = CSV.read("..\\data\\coalescence_data_slip_12_gamma0_1e-5_hmin_012.csv", DataFrame)
coalescene_data2 = CSV.read(HTTP.get("https://raw.githubusercontent.com/Zitzeronion/TimeDepWettabilityPaper/Drop_coalescence/Data_CSV/coalescence_data_slip_12_gamma0_1e-5_hmin_012.csv").body, DataFrame)

# ╔═╡ ee0a7ed7-9893-42f9-adfe-c23ff5315db0
md"To this dataframe we add the newly computed 
```math
	\Delta h = h_d^l - h_d^r
```
field.
"

# ╔═╡ fde19bfb-59ee-45bc-80e2-3c445ca18ce1
begin
	coalescene_data[!, "drop_diff"] .= abs.(coalescene_data.height_droplet_left .- coalescene_data.height_droplet_right)
	coalescene_data2[!, "drop_diff"] .= abs.(coalescene_data2.height_droplet_left .- coalescene_data2.height_droplet_right)
end

# ╔═╡ 350b5fa7-6b92-4f7b-b23f-a65b201c036e
md"Normalizing the lattice Boltzmann time steps we use the inertio capillary time scale.
It is computed according to 
```math
\tau_{ic} = \sqrt{\frac{\rho R^{3}}{\gamma}},
```
with $R$ being the initial radius and $\rho$ being the density, which we set 1."

# ╔═╡ d8f60516-34d7-4bd1-a330-d8a3becfc3e6
function tau_ic(;ρ=1, r₀=171, γ=1e-5)
	return sqrt(ρ*r₀^3/γ)
end

# ╔═╡ cc8ed16f-a3e2-4103-a03c-4414574f0921
md"Using this factor we create a new column called `t_norm` which contains the normalized time."

# ╔═╡ b3c4efff-2b13-4602-8df7-80d35af4af68
coalescene_data.t_norm .= coalescene_data.time ./ tau_ic() 

# ╔═╡ d6451808-f53c-4e98-82ad-30a79e697120
md"Below is the list of names for the different surface tension gardients.
On top a time set which looks reasonable spaced in loglog plots and does not collapse too much in the late simulation times.

And so more convenience dicts to store things." 

# ╔═╡ f184aa9f-2e04-4d22-b2fb-6f3297ea310e
begin
	# Surface tension names
	γ_names = ["const", "step", "tanh1", "tanh2", "tanh5","tanh10", "tanh20", "tanh50", "tanh100", "tanh200"]
	# Good subset of times for loglog plots
	log_t = [1, 2, 3, 4, 6, 10, 20, 30, 40, 60, 100, 200, 300, 400, 600, 900, 1001, 1002, 1003, 1005, 1007, 1010, 1015, 1022, 1030, 1040, 1055, 1078, 1110, 1160, 1250, 1360, 1500, 1700, 2000, 2500, 3200, 4200, 5600, 7800, 10900]
	# Labels for the different surface tension gradients
	label_dict = Dict("const" => "γ=γ₀", "step" => "γ=Θ(x)", "tanh1" =>"γ=s(x;1)", "tanh2" =>"γ=s(x;2)", "tanh5" =>"γ=s(x;5)", "tanh10" =>"γ=s(x;10)", "tanh20" =>"γ=s(x;20)", "tanh50" =>"γ=s(x;50)", "tanh100" =>"γ=s(x;100)", "tanh200" =>"γ=s(x;200)")
	# Some markers
	marker_dict = Dict(1 => :d, 2 => :s, 3 => :r, 4 => :h, 5 => :star4, 6 => :ut, 7 => :dt, 8 => :p, 9 => :lt, 10 => :s, 2 => :s, 2 => :s, 2 => :s)
end

# ╔═╡ b993564a-077d-47ed-b8b9-6548d2f0fbca
md"""
## Interactive

Now we can finally address the question on how $\Delta h$ changes based on the surface tension field.
In the plot below we just use three different lines, the constant field, the step and a smoothed variant.

The smoothing can be increased by moving the slider to the right

$(@bind l Slider(3:10))
"""

# ╔═╡ 06f3cb7a-544a-4065-85d0-7ffe76cf5105
begin
	i = 1
	plot(@subset(coalescene_data, :g_x .== γ_names[i]).t_norm, #[log_t] 
		@subset(coalescene_data2, :g_x .== γ_names[i]).drop_diff, 
				label=label_dict[γ_names[i]],
				# st = :scatter,
				l = (3, :solid),
				ylabel = "Δh", 
				xlabel = "t/τ",
				legend = :topleft,
				#yaxis = :log,
				grid = false,
				legendfontsize = 12,		# legend font size
    			tickfontsize = 14,			# tick font and size
    			guidefontsize = 15,
				# marker = (:circle, 8, 0.6, Plots.stroke(0, :gray)),
				)
	k = 2
	plot!(@subset(coalescene_data, :g_x .== γ_names[k]).t_norm, #[log_t] 
		@subset(coalescene_data, :g_x .== γ_names[k]).drop_diff, 
				label=label_dict[γ_names[k]],
				l = (3, :dash),
				# st = :scatter,
				# marker = (marker_dict[k-1], 8, 0.6, Plots.stroke(0, :gray)),
	)
	plot!(@subset(coalescene_data, :g_x .== γ_names[l]).t_norm, #[log_t] 
		@subset(coalescene_data, :g_x .== γ_names[l]).drop_diff, 
				label=label_dict[γ_names[l]],
				l = (3, :dashdot), 
				# st = :scatter,
				# marker = (marker_dict[l-1], 8, 0.6, Plots.stroke(0, :gray)),
	)
end

# ╔═╡ 8c713445-d87a-47e2-a753-e6f6f4812a7a
md"Having tested some smoothing widths and learned how they affect the droplets size we like to collect this results in a plot.
In the plot we use a subset of smoothing widths and plot the evolution of the height difference.
"

# ╔═╡ 8ffb2721-e939-45b0-8566-ad8c350a05b4
@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 10)
    n = length(y)
    sx, sy = x[1:step:n], y[1:step:n]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := [Inf]
        y := [Inf]
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end

# ╔═╡ 8707f1b4-8418-4ebc-99e1-b6aba22c25d2
begin
	colors_h = [palette(:default)[5], palette(:default)[6], palette(:default)[7], palette(:default)[8],palette(:default)[9]]
	hdiff = plot(@subset(coalescene_data, :g_x .== γ_names[i]).t_norm, #[log_t] 
		@subset(coalescene_data2, :g_x .== γ_names[i]).drop_diff, 
				label=label_dict[γ_names[i]],
				# st = :scatter,
				l = (3, :solid),
				ylabel = "Δh", 
				xlabel = "t/τ",
				legend = :topleft,
				# yaxis = :log,
				grid = false,				
		      	st = :samplemarkers, 				# some recipy stuff
		        step = 1000, 						# density of markers
		        marker = (:circle, 8, 0.6, Plots.stroke(0, :gray)),	
				legendfontsize = 12,		# legend font size
    			tickfontsize = 14,			# tick font and size
    			guidefontsize = 15,
				# ,
				)
	plot!(@subset(coalescene_data, :g_x .== γ_names[2]).t_norm, #[log_t] 
			@subset(coalescene_data, :g_x .== γ_names[2]).drop_diff, 
					label=label_dict[γ_names[2]],
					l = (3, :dash),
					st = :samplemarkers, 				# some recipy stuff
		        	step = 1000,
					marker = (marker_dict[k-1], 8, 0.6, Plots.stroke(0, :gray)),
		)
	for k in 5:7
		plot!(@subset(coalescene_data, :g_x .== γ_names[k]).t_norm, #[log_t] 
			@subset(coalescene_data, :g_x .== γ_names[k]).drop_diff, 
					label=label_dict[γ_names[k]],
					l = (3, :auto, colors_h[k-4]),
					st = :samplemarkers, 				# some recipy stuff
		        	step = 1000,
					# st = :scatter,
					marker = (marker_dict[k-1], colors_h[k-4], 8, 0.6, Plots.stroke(0, :gray)),
		)
	end
	plot!(@subset(coalescene_data, :g_x .== γ_names[8]).t_norm, #[log_t] 
			@subset(coalescene_data, :g_x .== γ_names[8]).drop_diff, 
					label=label_dict[γ_names[8]],
					l = (3, :auto, colors_h[4]),
					st = :samplemarkers, 				# some recipy stuff
		        	step = 1000,
					marker = (:dt, colors_h[4], 8, 0.6, Plots.stroke(0, :gray)),
					# ylims=(1e-3, 100)
		)
end

# ╔═╡ e2dd8091-3a5e-448c-87df-f8bdc734afdc
begin
	plot(@subset(coalescene_data, :g_x .== γ_names[i]).t_norm, #[log_t] 
		@subset(coalescene_data2, :g_x .== γ_names[i]).drop_diff, 
				label=label_dict[γ_names[i]],
				# st = :scatter,
				l = (3, :solid),
				ylabel = "Δh", 
				xlabel = "t/τ",
				legend = :bottomright,
				axis = :log,
				minorticks=10,
				grid = false,				
		      	st = :samplemarkers, 				# some recipy stuff
		        step = 1000, 						# density of markers
		        marker = (:circle, 8, 0.6, Plots.stroke(0, :gray)),	
				legendfontsize = 12,		# legend font size
    			tickfontsize = 14,			# tick font and size
    			guidefontsize = 15,
				xlims= (10, 200),
				ylims = (0.1, 20)
				)
	plot!(@subset(coalescene_data, :g_x .== γ_names[2]).t_norm, #[log_t] 
			@subset(coalescene_data, :g_x .== γ_names[2]).drop_diff, 
					label=label_dict[γ_names[2]],
					l = (3, :dash),
					st = :samplemarkers, 				# some recipy stuff
		        	step = 1000,
					marker = (marker_dict[k-1], 8, 0.6, Plots.stroke(0, :gray)),
		)
	for k in 5:7
		plot!(@subset(coalescene_data, :g_x .== γ_names[k]).t_norm, #[log_t] 
			@subset(coalescene_data, :g_x .== γ_names[k]).drop_diff, 
					label=label_dict[γ_names[k]],
					l = (3, :auto),
					st = :samplemarkers, 				# some recipy stuff
		        	step = 1000,
					# st = :scatter,
					marker = (marker_dict[k-1], 8, 0.6, Plots.stroke(0, :gray)),
		)
	end
	plot!(@subset(coalescene_data, :g_x .== γ_names[8]).t_norm, #[log_t] 
			@subset(coalescene_data, :g_x .== γ_names[8]).drop_diff, 
					label=label_dict[γ_names[8]],
					l = (3, :auto),
					st = :samplemarkers, 				# some recipy stuff
		        	step = 1000,
					marker = (marker_dict[8], 8, 0.6, Plots.stroke(0, :gray)),
					# ylims=(1e-3, 100)
		)
	plot!(collect(1:0.1:300), 0.2 * collect(1:0.1:300).^(3/4))
end

# ╔═╡ 8377837a-24b9-47ad-9777-1f9c89c32d41
# pwd()
savefig(hdiff, "..\\figures\\hdiff.svg")

# ╔═╡ 42dc5c99-1980-421e-b993-2ea92e6fde5c
md"### Separation time

One further question we can ask the data is when the droplets separate.
The separation can be defined using
```math
	\min_t(h_0) \leq h_{\ast},
```
meaning when $h_0(t) \leq h_{\ast}$ we assume the droplets are separated by an *dry* spot. 

This only happens for a few of our simulations, in fact only for $w<50$.
What we know so far is when the dip of the bridge starts to build up.
But that dip is correlated with the separation time.
However because set a maximal simulation time, we do not observe every separation. 
On the other hand, for large values of $w > 50$ it is more likely that the other droplet is just eaten up.
That is why the plot only shows a subset of smearing widths.
"

# ╔═╡ 0c175894-a5dd-4844-af32-e24a2d1dc03a
begin
	sep_time = Float64[]
	gamma_n = Int[]
	for i in γ_names[3:end]
		one_gamma = @subset(coalescene_data, :g_x .== i)
		if minimum(one_gamma.bridge_min) < 0.09
			# println("Hello, $i")
			push!(sep_time, first(one_gamma[one_gamma.bridge_min .≤ 0.09, :t_norm]))
		else
			push!(sep_time, 200)
		end
		# println(first(one_gamma[one_gamma.bridge_min .≤ 0.09, :t_norm]))
		# println(i[5:end])
		push!(gamma_n, parse(Int, i[5:end]))
	end
end

# ╔═╡ a2ed13b5-84db-4f02-bd6d-39f39ba5642c
begin
	tau_sep = plot(gamma_n[1:5] ./ 171, sep_time[1:5],
		st = :scatter,
		xlabel = "w/R₀",
		ylabel = "t/τ",
		grid=:false,
		label="",
		marker = (:circle, 30, 0.6, Plots.stroke(0, :gray)),
		xticks = ([0.00, 0.04, 0.08, 0.12], ["0.00", "0.04", "0.08", "0.12"]),
		legendfontsize = 20,		# legend font size
    	tickfontsize = 18,			# tick font and size
    	guidefontsize = 20,
		legend=:topleft,
		xlims=(0,0.13),
		ylims=(0,80),
		)
	data_x = collect(0:0.001:1)
	Δγ = 8e-6
	plot!(data_x , data_x.^(3/2) .* 1700, l=(5, :dash, :black), label="")
	# plot!(data_x , data_x.^2 ./ (22*Δγ), l=(3, :dash, :gray), label="w²/(2(γ₀-Δγ)")
end

# ╔═╡ ba644615-6717-4598-81cc-b8484e5f4855
savefig(tau_sep, "..\\figures\\seperation_time_scaling.svg")

# ╔═╡ e6a1d40d-3ae1-4608-9788-e6a149d8f283
md"What is interesting in this plot is the linear depenedency on the smearing width.

We can use some function
```math
	y(w) = k\cdot w + w_0
```
where $k$ is a factor for the slope and $x_0$ is an offset.
Clearly the offset is small, so $w_0 = 0$ for now.
More interesting is the slope.
For now I have no understanding why the slope takes the value it has.
Problem for Stefan after the reviews.
"

# ╔═╡ 7d4c9009-ac77-4a4d-9010-5c69681b2004
md"## Bridge velocity

In our simulations we saw that the bridge is moving.
Meaning that changes it's position with time.
While Marangoni flow drives fluid into regions of higher surface tension we observe a motion of the bridge downstream.
Fluid passes into the left droplet and increases the bridge height $h_0$, by doing so it pushes the bridge to the right.
Until the bridge does not longer exists and the droplets are either separated or coalesced.

One hypothesis is that the droplets separate when the bridge traveled to a region without surface tension gradient, 

```math
	x_{sep} := \gamma(x) = \gamma_0 - \Delta\gamma .
```

This can be tested easily as we have data on the bridge position for all runs.
Let's plot it for a few examples
"

# ╔═╡ 7eff52dc-5ae5-4655-b8ad-1e366b61dbbc
begin
	#plot(@subset(coalescene_data, :g_x .== γ_names[i]).t_norm, #[log_t] 
		# @subset(coalescene_data2, :g_x .== γ_names[i]).neck_pos .- 172, 
				# label=label_dict[γ_names[i]],
				# st = :scatter,
				# l = (3, :solid),
				# ylabel = "χ/R₀", 
				# xlabel = "t/τ",
				# legend = :topleft,
				# xaxis = :log,
				# grid = false,				
		      	# st = :samplemarkers, 				# some recipy stuff
		        # step = 1000, 						# density of markers
		        # marker = (:circle, 8, 0.6, Plots.stroke(0, :gray)),	
				# legendfontsize = 12,		# legend font size
    			# tickfontsize = 14,			# tick font and size
    			# guidefontsize = 15,
				# xlim = (0.1, 100),
				# ylim = (-2, 30),
				# )
	bv = plot(@subset(coalescene_data, :g_x .== γ_names[3]).t_norm, #[log_t] 
		(@subset(coalescene_data, :g_x .== γ_names[3]).neck_pos .- 512) ./ 171, 
				label=label_dict[γ_names[3]],
				# st = :scatter,
				l = (3, :solid),
				ylabel = "χ/R₀", 
				xlabel = "t/τ",
				legend = :topleft,
				axis = :log,
				grid = false,				
		      	# st = :samplemarkers, 				# some recipy stuff
		        # step = 1000, 						# density of markers
		        # marker = (:circle, 8, 0.6, Plots.stroke(0, :gray)),	
				legendfontsize = 12,		# legend font size
    			tickfontsize = 14,			# tick font and size
    			guidefontsize = 15,
				xlim = (0.1, 100),
				# ylim = (-2, 30),
				)

	plot!(@subset(coalescene_data, :g_x .== γ_names[5]).t_norm, #[log_t] 
		(@subset(coalescene_data, :g_x .== γ_names[5]).neck_pos .- 512) ./ 171, 
				label=label_dict[γ_names[5]],
				# st = :scatter,
				l = (3, :solid),
				)
	plot!(@subset(coalescene_data, :g_x .== γ_names[6]).t_norm, #[log_t] 
		(@subset(coalescene_data, :g_x .== γ_names[6]).neck_pos .- 512) ./ 171, 
				label=label_dict[γ_names[6]],
				# st = :scatter,
				l = (3, :solid),
				)
	plot!(@subset(coalescene_data, :g_x .== γ_names[7]).t_norm, #[log_t] 
		(@subset(coalescene_data, :g_x .== γ_names[7]).neck_pos .-512) ./ 171, 
				label=label_dict[γ_names[7]],
				# st = :scatter,
				l = (3, :solid),
				)
	# scatter!(sep_time[3], 512 .- @subset(coalescene_data, :g_x .== γ_names[5]).neck_pos[])
end

# ╔═╡ b428d804-6211-45be-b753-41a7580f1d52
savefig(bv, "..\\figures\\bridge_vel_log.svg")

# ╔═╡ 7ffa2d9f-8bee-4775-9f65-81fe05875e1d
begin
	for i in 1:5
		ts = i
		nm = ts+2
		println(coalescene_data[(coalescene_data.g_x .== γ_names[nm]) .& (coalescene_data.t_norm .== sep_time[ts]), :neck_pos])
	end
end

# ╔═╡ c02d1b7a-b8f8-457c-96f2-383a6422a36d
md"### Pressures

One last thing that we are looking into is the distribution of the pressures.
The pressure is something similar to 

```math
	p = \gamma(\Delta h - \Pi(h)),
```
therefore a second derivative of the height with an additional disjoining pressure component.
One natural unit of pressure for the coalescence is

```math
	p_c = \frac{\gamma}{R},
```
often called capillary pressure.

In the following we take a look at a subset of simulations where we measured the pressures and compare them to the capillary pressure.
"

# ╔═╡ ff16591c-26a8-4330-96d8-dd807b5bc9fb
pressure_data = CSV.read("..\\data\\pressure_data.csv", DataFrame)

# ╔═╡ 7a244fc8-a74d-4e12-892e-570cef474bca
begin
	xaxis = -171:2:171 
	pressure_plot = plot()
	# for i in [10000, 1000000, 100000000]
	ci = [palette(:default)[1], palette(:default)[2], palette(:default)[3]] 
	for i in enumerate(["step", "tanh5", "tanh100"])
		plot!(xaxis ./171, 
			pressure_data[(pressure_data.kind .== "lap") .& (pressure_data.time .== 100000000) .& (pressure_data.gamma .== i[2]), :pressure],
			l=(3, :solid, ci[i[1]]),
			label="$(label_dict[i[2]])",
			xlabel="x/R₀",
			ylabel="P/p₀",
			st = :samplemarkers, 				# some recipy stuff
			step = 11, 
			marker = (:circle, 8, 0.6, Plots.stroke(0, :gray), ci[i[1]]),
			legendfontsize = 12,		# legend font size
	    	tickfontsize = 14,			# tick font and size
	    	guidefontsize = 15,
			legend=:bottomleft,
			grid=false,
		)
		plot!(xaxis ./171, 
			pressure_data[(pressure_data.kind .== "disj") .& (pressure_data.time .== 100000000) .& (pressure_data.gamma .== i[2]), :pressure],
			l=(3, :dash, ci[i[1]]),
			label="",
			st = :samplemarkers, 				# some recipy stuff
			step = 11, 
			marker = (:star, 8, 0.6, Plots.stroke(0, :gray), ci[i[1]]),
		)
	end
	# lens!([-0.15,0.85], [-1.5, 1.5], inset = (1, bbox(0.1, 0.2, 0.4, 0.4)))
	plot!(xlims=(-0.85, 0.65))
	lens!([-0.15,0.85], [-0.5, 0.5], grid=false, inset = (1, bbox(0.1, 0.2, 0.4, 0.4)))
end

# ╔═╡ a30ea74e-3614-4eb8-91f4-5d89a0f2f73c
savefig(pressure_plot, "..\\Figures\\pressures_diff_gam.svg")

# ╔═╡ d3e3929c-97e5-4635-ad25-54ea87c24c82
md"## Interfaces

With all the data it is actually helpful to take a look at the fluid-vapor interface.
We know there is about three different scenarios
- Coalescence
- Separation
- Asymmetric coalescence
" 

# ╔═╡ 4eb19849-a429-4a69-a382-f4e04e706105
begin
	tt = tau_ic()
	go_back = "C:\\Users\\zitz\\Software_Projects\\Swalbe.jl\\data\\Drop_coalescence_long"
    folder = "\\gamma_10_"
	times_plot = []
	surfs = ["const", "tanh5", "tanh100"]
	snapshots = zeros(1024, 3, 3)
	for i in enumerate(surfs)
		for t in enumerate([10000, 10000000, 100000000])
			push!(times_plot, t[2] / tt)
			df = load(go_back * folder * "$(i[2])_periodic_tmax_100000000_slip_12_L_1024_hm_12_hc_3_gamma_10.jld2") |> DataFrame
		snapshots[:, i[1], t[1]] .= df[!, Symbol("h_$(t[2])")]
		end
	end
end

# ╔═╡ 6e01e3cf-1f51-4cff-a87c-0c8c524d415a
begin
	lfs = 14
	tfs = 18
	gfs = 20
	p1 = plot(collect(-511:512) ./ 171, snapshots[:, 1, 1] ./ 171, 
		annotations = (1.8, 0.23, Plots.text("(a)", 24, :left)),
		label="t = $(round(times_plot[1],sigdigits=2))τ",
		xlabel="x/R₀", ylabel="h/R₀",
		st = :samplemarkers,
		w=3,
		c=ci[1],
		step = 40, 
		marker = (:circle, 8, 0.6, Plots.stroke(0, :gray), ci[1]),
		legendfontsize = lfs,		# legend font size
	    tickfontsize = tfs,			# tick font and size
	    guidefontsize = gfs,
		legend=:topleft,
		grid=false,
		)
	plot!(collect(-511:512) ./ 171, snapshots[:, 1, 2] ./ 171, 
		label="t = $(round(times_plot[2],sigdigits=2))τ",
		st = :samplemarkers,
		l=(3, :solid, ci[1]),
		step = 40, 
		marker = (:star, 8, 0.6, Plots.stroke(0, :gray), ci[1]),
		)
	plot!(collect(-511:512) ./ 171, snapshots[:, 1, 3] ./ 171, 
		label="t = $(round(times_plot[3],sigdigits=2))τ",
		st = :samplemarkers,
		l=(3, :solid, ci[1]),
		step = 40, 
		marker = (:ut, 8, 0.6, Plots.stroke(0, :gray), ci[1]),
		)
	ylims!(0, 0.25)
	xlims!(-2.1, 2.1)
end

# ╔═╡ 02fad95f-1e0f-4444-9a31-6b6aedd299fe
begin
	p2 = plot(collect(-511:512) ./ 171, snapshots[:, 2, 1] ./ 171,xlabel="x/R₀", 
			label="",
		annotations = (1.8, 0.23, Plots.text("(b)", 24, :left)),
			st=:samplemarkers,
			l=(3, ci[2]),
			step = 40, 
		legendfontsize = lfs,		# legend font size
	    tickfontsize = tfs,			# tick font and size
	    guidefontsize = gfs,
		legend=:topleft,
		grid=false,
			marker = (:circle, 8, 0.6, Plots.stroke(0, :gray), ci[2]))
		plot!(collect(-511:512) ./ 171, snapshots[:, 2, 2] ./ 171, 
			label="",
			st=:samplemarkers,
			l=(3, ci[2]),
			step = 40, 
			marker = (:star, 8, 0.6, Plots.stroke(0, :gray), ci[2]))
		plot!(collect(-511:512) ./ 171, snapshots[:, 2, 3] ./ 171, 
			label="",
			st=:samplemarkers,
			l=(3, ci[2]),
			step = 40, 
			marker = (:ut, 8, 0.6, Plots.stroke(0, :gray), ci[2]))
	ylims!(0, 0.25)
	xlims!(-2.1, 2.1)
end

# ╔═╡ b8fd374a-ec9d-4f65-97e4-d02c77ff57ea
begin
	p3 = plot(collect(-511:512) ./ 171, snapshots[:, 3, 1] ./ 171, xlabel="x/R₀",
		label="",
		annotations = (1.8, 0.23, Plots.text("(c)", 24, :left)),
			st=:samplemarkers,
			l=(3, ci[3]),
			step = 40, 
		legendfontsize = lfs,		# legend font size
	    tickfontsize = tfs,			# tick font and size
	    guidefontsize = gfs,
		legend=:topleft,
		grid=false,
			marker = (:circle, 8, 0.6, Plots.stroke(0, :gray), ci[3]))
		plot!(collect(-511:512) ./ 171, snapshots[:, 3, 2] ./ 171, label="",
			st=:samplemarkers,
			l=(3, ci[3]),
			step = 40, 
			marker = (:star, 8, 0.6, Plots.stroke(0, :gray), ci[3]))
		plot!(collect(-511:512) ./ 171, snapshots[:, 3, 3] ./ 171, label="",
			st=:samplemarkers,
			l=(3, ci[3]),
			step = 40, 
			marker = (:ut, 8, 0.6, Plots.stroke(0, :gray), ci[3]))
	ylims!(0, 0.25)
	xlims!(-2.1, 2.1)
end

# ╔═╡ 663f5e86-9008-425f-a1b8-41d0231b1576
begin
	savefig(p1, "..\\Figures\\h_final_three_1.svg")
	savefig(p2, "..\\Figures\\h_final_three_2.svg")
	savefig(p3, "..\\Figures\\h_final_three_3.svg")
end

# ╔═╡ 57f46a0f-080b-4d14-bdb0-329f04f424c2
md"## Surface tension fields

For a better understanding of the surface tension fields we supply a plot of the various $\gamma(x)$.

We use three different fields for this study.
1. Constant surface tension
2. A surface tension step
3. Smoothing using a $\tanh$ with variable smearing widths $w$

"

# ╔═╡ 4e7a81c7-e8ac-4733-862c-8de6d7589faa
begin 
	L = 1024
	γ₀ = 1e-5
end

# ╔═╡ 20296440-d014-4ce3-9f56-ab91393f9b9e
"""
	const_gamma()

Constant surface tension over the whole domain.
"""
function const_gamma(;L=L, γ=γ₀)
    return fill(γ, L)
end

# ╔═╡ 2321e432-65be-4df6-80f2-c91eef4e9a84
"""
	step_gamma()

A step like surface tension field, using the equation below

`` \\gamma(x) = \\begin{cases} \\gamma_0\\qquad\\qquad\\text{for}~x < L/2 \\\\ \\gamma_0 - \\frac{2 \\gamma_0}{10}\\quad\\text{else}
	\\end{cases}
``
"""
function step_gamma(;L=L, γ=γ₀, perc=20)
    x = ones(L)
    for i in 1:L
        if i < L÷2
            x[i] = γ
        else
            x[i] = γ - (γ * perc / 100)
        end
    end
    return x
end

# ╔═╡ 4a287282-6aae-491e-a8f9-40351b7d8d7c
"""
	tanh_gamma()

A tangent hyperbolicus to smooth the surface tension step, using the equation

`` \\gamma(x) = \\gamma_0\\cdot \\bigg|1 - \\big( \\frac{1}{2} - s(x,l,w) \\big)\\bigg| + ( \\gamma_0 - \\Delta \\gamma) \\bigg[1 - \\bigg|1 - \\big( \\frac{1}{2} - s(x,l,w) \\big) \\bigg| \\bigg]   
``
"""
function tanh_gamma(;L=L, γ=γ₀, perc=20, sl=1)
    l = collect(1.0:L)
    function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	x = ones(L)
	x .= γ .* smooth(l, L, sl) .+ (1 .- smooth(l, L, sl)) .* (γ - (γ * perc / 100)) 
	return x
end

# ╔═╡ 22a15d69-4435-48d1-90a8-9187987e49da
begin
	x_r0 = collect(-511:512) ./ 171
	# step_gamma(γ=s), tanh_gamma(sl=1, γ=s),
	Gammas = plot(x_r0, const_gamma(γ=γ₀) ./ const_gamma(γ=γ₀),  label="γ(x) = γ₀",
		xlabel="x/R₀", ylabel="γ/γ₀",
		st = :samplemarkers,
		l=(3),
		step = 40, 
		marker = (:circle, 8, 0.6, Plots.stroke(0, :gray)),
		legendfontsize = 12,		# legend font size
	    tickfontsize = 14,			# tick font and size
	    guidefontsize = 15,
		legend=:bottomleft,
		grid=false,
		ylim=(0.78, 1.02))
	plot!(x_r0, step_gamma(γ=γ₀) ./ const_gamma(γ=γ₀),  label="γ(x) = Θ(x)",
		st = :samplemarkers,
		l=(3, :dash),
		step = 40, 
		marker = (:d, 8, 0.6, Plots.stroke(0, :gray)))
	plot!(x_r0, tanh_gamma(γ=γ₀, sl=5) ./ const_gamma(γ=γ₀),  label="γ(x) = s(x;5)",
		st = :samplemarkers,
		step = 40, 
		l=(3, palette(:default)[5], :dashdot),
		marker = (:hex, 8, palette(:default)[5], 0.6, Plots.stroke(0, :gray)))
	plot!(x_r0, tanh_gamma(γ=γ₀, sl=100) ./ const_gamma(γ=γ₀),  label="γ(x) = s(x;100)",
		st = :samplemarkers,
		l=(3, palette(:default)[9], :dot),
		step = 40, 
		marker = (:p, 8, 0.6, palette(:default)[9], Plots.stroke(0, :gray)))
end

# ╔═╡ bd2b3f2f-cfdf-4ff4-ba59-0353a91f07ea
savefig(Gammas, "..\\Figures\\gammas.svg")

# ╔═╡ 3cd97467-9991-485b-8e33-a07655038c8e
begin
	x = tanh_gamma(γ=γ₀, sl=5)
	plot(x, xlim = (500, 524))
	grads = zeros(1024)
	grads .= circshift(x, 1) .- circshift(x, -1)
	grads[1:3] .= grads[1022:end] .= 0
	plot!(20 .* grads)
end

# ╔═╡ edde35ab-d4fe-4fc2-b63f-c1d198df6b99
md"
### Other Fluids

After this scientific journy into the world of droplets and their attraction let's end this notebook with something tasty.
There is this delicious mixture, called espresso tonic.
Interested?

Check out the video from James Hoffmann while I heat up the espresso machine.
"

# ╔═╡ 1429de49-123c-4ce9-8713-b654379f61f5
html"""

<div style="padding:10% 0 0 0;position:relative;"><iframe width="560" height="315" src="https://www.youtube.com/embed/4x92J4gQwzM" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
DataFramesMeta = "~0.11.0"
FileIO = "~1.14.0"
HTTP = "~0.9.17"
JLD2 = "~0.4.22"
Plots = "~1.29.1"
PlutoUI = "~0.7.39"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "7297381ccb5df764549818d9a7d57e45f1057d30"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.18.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "0f4e115f6f34bbe43c19751c90a38b2f380637b9"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.3"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "f1d89a07475dc4b03c08543d1c6b4b2945f33eca"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.11.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c98aea696662d09e215ef7cda5296024a9646c75"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "3a233eeeb2ca45842fe100e0413936834215abf5"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.4+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "9e42de869561d6bdf8602c57ec557d43538a92f0"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.29.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2bbd9f2e40afd197a1379aef05e0d85dba649951"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2c11d7290036fe7aac9038ff312d3b3a2a5bf89e"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.4.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "9abba8f8fb8458e9adf07c8a2377a070674a24f1"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.8"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─7ed11bdd-2f22-44a5-94c2-8ef0d841769e
# ╠═50901b90-eb08-11ec-3087-81a949a22d5d
# ╟─cb1d2d4f-5f05-4fe9-93db-973c9113cf67
# ╠═ef661fd9-e3c2-4fd3-a36b-75a43222155f
# ╠═0671ef84-5f84-4967-bf29-b3f9e83f9f08
# ╟─ee0a7ed7-9893-42f9-adfe-c23ff5315db0
# ╟─fde19bfb-59ee-45bc-80e2-3c445ca18ce1
# ╟─350b5fa7-6b92-4f7b-b23f-a65b201c036e
# ╠═d8f60516-34d7-4bd1-a330-d8a3becfc3e6
# ╟─cc8ed16f-a3e2-4103-a03c-4414574f0921
# ╠═b3c4efff-2b13-4602-8df7-80d35af4af68
# ╟─d6451808-f53c-4e98-82ad-30a79e697120
# ╠═f184aa9f-2e04-4d22-b2fb-6f3297ea310e
# ╟─b993564a-077d-47ed-b8b9-6548d2f0fbca
# ╟─06f3cb7a-544a-4065-85d0-7ffe76cf5105
# ╟─8c713445-d87a-47e2-a753-e6f6f4812a7a
# ╟─8ffb2721-e939-45b0-8566-ad8c350a05b4
# ╠═8707f1b4-8418-4ebc-99e1-b6aba22c25d2
# ╠═e2dd8091-3a5e-448c-87df-f8bdc734afdc
# ╠═8377837a-24b9-47ad-9777-1f9c89c32d41
# ╟─42dc5c99-1980-421e-b993-2ea92e6fde5c
# ╠═0c175894-a5dd-4844-af32-e24a2d1dc03a
# ╠═a2ed13b5-84db-4f02-bd6d-39f39ba5642c
# ╠═ba644615-6717-4598-81cc-b8484e5f4855
# ╟─e6a1d40d-3ae1-4608-9788-e6a149d8f283
# ╠═7d4c9009-ac77-4a4d-9010-5c69681b2004
# ╠═7eff52dc-5ae5-4655-b8ad-1e366b61dbbc
# ╠═b428d804-6211-45be-b753-41a7580f1d52
# ╠═7ffa2d9f-8bee-4775-9f65-81fe05875e1d
# ╟─c02d1b7a-b8f8-457c-96f2-383a6422a36d
# ╠═ff16591c-26a8-4330-96d8-dd807b5bc9fb
# ╠═7a244fc8-a74d-4e12-892e-570cef474bca
# ╠═a30ea74e-3614-4eb8-91f4-5d89a0f2f73c
# ╟─d3e3929c-97e5-4635-ad25-54ea87c24c82
# ╠═4eb19849-a429-4a69-a382-f4e04e706105
# ╠═6e01e3cf-1f51-4cff-a87c-0c8c524d415a
# ╠═02fad95f-1e0f-4444-9a31-6b6aedd299fe
# ╠═b8fd374a-ec9d-4f65-97e4-d02c77ff57ea
# ╠═663f5e86-9008-425f-a1b8-41d0231b1576
# ╟─57f46a0f-080b-4d14-bdb0-329f04f424c2
# ╠═4e7a81c7-e8ac-4733-862c-8de6d7589faa
# ╟─20296440-d014-4ce3-9f56-ab91393f9b9e
# ╟─2321e432-65be-4df6-80f2-c91eef4e9a84
# ╟─4a287282-6aae-491e-a8f9-40351b7d8d7c
# ╠═22a15d69-4435-48d1-90a8-9187987e49da
# ╠═bd2b3f2f-cfdf-4ff4-ba59-0353a91f07ea
# ╠═3cd97467-9991-485b-8e33-a07655038c8e
# ╠═edde35ab-d4fe-4fc2-b63f-c1d198df6b99
# ╟─1429de49-123c-4ce9-8713-b654379f61f5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
