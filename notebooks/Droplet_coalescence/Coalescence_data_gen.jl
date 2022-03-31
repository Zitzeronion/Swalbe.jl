### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ a539f57c-7d60-4851-ba5e-60dad313fab7
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	
    using Swalbe
end

# ╔═╡ daa88477-ee4d-4563-a0ee-d228a672fd34
using Plots, Revise, DataFrames, FileIO, DataFramesMeta, StatsBase, CSV

# ╔═╡ 26520cda-a0ab-403c-9bbe-d39c93d85f7c
include("helpers.jl")

# ╔═╡ c18ba690-9f8e-11ec-1a41-7330ad3642ec
md"# Droplet coalescence - Data generation

This notebook is intended to create all data that we use to understand the coalescence of two sessile droplets on a surface tension gradient.
It uses `Swalbe.jl` to run the experiments or if the data is readly available skips the computation.

Depending on the simulation time the data generation can take a few hours.
**Be careful with the data** so that no compuation time is wasted.

In the following will be the creation of three data sets,

- A sample set
- Different surface tension gradients γ(x)
- Long runs with γ(x) ∝ tanh(f(x))
"

# ╔═╡ 499e69d3-f894-4c41-b320-021213184c3c
md"## First run

In this first run we asses that Swalbe.jl can perform simulation of noncoalescence.
We therefore perform up to 5×10⁶ time step with different surface tension gradients, 

- Constant γ
- Linear γ
- Step in γ
- Smoothed step with tanh of γ"

# ╔═╡ cdc87426-d742-4bf0-9699-c793e5753961
begin
	sys = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=5000000, δ=5.0))
	# Memory
	data = zeros(4, 50000, 1024)
	# Time
	sim_time = 100:100:size(data)[2]*100
	# Loop over γ
	for i in 1:4
		# Check if there is already a file created
		sim_name = "gamma_$(gamma_labels[i])_tmax_$(sys.param.Tmax).jld2"
		save_file = string(data_path, sim_name)
		# If so just read it from disc
		if isfile(save_file)
			println("There is already a file at location $(save_file)\nIf you still want to run the experiment change `sim_name`, can take several minutes.")
		# If not, compute the evolution of the droplet coalescence
		else
	 		data[i, :, :] = Swalbe.run_gamma(sys, γ[i, :], r₁=rad, r₂=rad)
			# Think about storing them on the disc
			df_fluid = Dict()
			# Loop through the time once more
			for t in 1:size(data)[2]
			    df_fluid["h_$(sim_time[t])"] = data[i,t,:]
			end
			# Save it
			println("Writing file $(gamma_labels[i]) to disk")
			save(save_file, df_fluid)
		end
		# Print that you are done
		println("Done with iteration $(gamma_labels[i])")
	end
end

# ╔═╡ 677a3ed5-046e-49de-a2bf-7f65fe904be4
md"## Second run

With the information the first run supplied we know that droplets stay seperated.
We now like to understand the noncoalescence a little bit better.
For that reason we focus on the smoothed step with tanh and take a look on the dependence of the smoothing width.
"

# ╔═╡ 0f007740-d067-4b89-9aeb-2b63f7c77d09
md"## Third run

The smoothing does in fact some interesting things.
Independent of the width the bridge starts to grow.
However depending on the width there is a time where the drops suddenly seperate.
We want to study this in more detail and for even longer time scales.

My guess there is that there should be a critical bridge height where the two can not sperate anymore, but this is just a guess.
"

# ╔═╡ f6ebd1b2-fbd3-4622-ab9d-dc98f50bb31f
begin
	sys_tanh_l = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=50000000, δ=5.0, tdump=1000))
	tanh_l_dict = Dict(1 => "sw_20", 2 => "sw_30", 3 => "sw_40", 4 => "sw_50")
	tanh_v_dict = Dict(1 => 20, 2 => 30, 3 => 40, 4 => 50)
	data_t = zeros(length(tanh_l_dict), sys_tanh_l.param.Tmax÷sys_tanh_l.param.tdump, sys_tanh_l.L)
	# Memory
	nexp = length(tanh_l_dict)
	data_tanh_l = zeros(nexp, sys_tanh_l.param.Tmax÷sys_tanh_l.param.tdump, sys_tanh_l.L)
	# Time here
	sim_time_t = 1000:1000:50000000
	# Loop over γ
	for i in 1:nexp
		ggs = zeros(1024)
		gamma_curves_tanh!(ggs, sl=tanh_v_dict[i])
		# Check if there is already a file created
		sim_name_tanh = "gamma_tanh_width_$(tanh_l_dict[i])_tmax_$(sys_tanh_l.param.Tmax).jld2"
		save_file = string(data_path, sim_name_tanh)
		# If so just read it from disc
		if isfile(save_file)
			println("There is already a file at location $(save_file)\nIf you still want to run the experiment change `sim_name_tanh`, can take several minutes.")
			
		# If not, compute the evolution of the droplet coalescence
		else
	 		data_t[i, :, :] = Swalbe.run_gamma(sys_tanh_l, ggs, r₁=rad, r₂=rad, dump=sys_tanh_l.param.tdump)
			# Think about storing them on the disc
			df_fluid = Dict()
			# Loop through the time once more
			for t in 1:size(data_t)[2]
		        df_fluid["h_$(sim_time_t[t])"] = data_t[i,t,:]
		    end
			# Save it
			println("Writing file $(tanh_l_dict[i]) to disk")
			save(save_file, df_fluid)
		end
		# Print that you are done
		println("Done with iteration $(tanh_l_dict[i])")
	end
end

# ╔═╡ 377a8f41-fc08-4941-8b08-530da4f0a987
md"## Forth run

Since there is an effect due to the periodicity of the domain we like to experiment a little with walls.
For this reason we use bounce back boundaries at the start and end of the domain.
This will probably lead to some minuscus at the boundary but nothing to bad happens hopefully.
"

# ╔═╡ 604ce508-e7a6-422f-8a55-e1658143bbac
function run_gamma_walls(
    sys::Swalbe.SysConstWithBound_1D,
    gamma::Vector;
    r₁=115,
    r₂=115, 
    θ₀=1/9,  
    verbos=true, 
    dump = 100, 
    fluid=zeros(sys.param.Tmax÷dump, sys.L)
)
    println("Simulating droplet coalecense with surface tension gardient on a domain with walls")
	Swalbe.obslist!(sys)
    state = Swalbe.Sys(sys, kind="gamma_bound")
    drop_cent = (sys.L/3, 2*sys.L/3)
    state.basestate.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state, sys)
    state.γ .= gamma
	Swalbe.∇γ!(state)
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.param.Tmax
        if verbos
            if t % sys.param.tdump == 0
                mass = 0.0
                mass = sum(state.basestate.height)
                println("Time step $t bridge height is $(round(minimum(state.basestate.height[sys.L÷2-50:sys.L÷2+50]), digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, γ=gamma)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.basestate.F .= -state.basestate.h∇p .- state.basestate.slip .- state.∇γ
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(fluid, state.basestate.height, t, dumping = dump)
    end
	
    return fluid
    
end

# ╔═╡ d7e92ccc-a1c5-4547-a6f6-5311514962a1
begin
	L = 1024
	obst = zeros(L)
    obst[1]=1
    obst[2]=1
    obst[L]=1
    obst[L-1]=1
    sys_bound = Swalbe.SysConstWithBound_1D{Float64}(obs=obst, L=L,param= Swalbe.Taumucs(Tmax=1000000, n=9, m=3))
	gamma = 0.0001 .* (1 .- 0.1 .* collect(1.0:L) ./ L)
	fl = run_gamma_walls(sys_bound, gamma, r₁=500, r₂=500)
end

# ╔═╡ e073abd0-dcd2-4c3d-a043-238f81952446
begin
	sys_tanh = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=10000000, δ=5.0, tdump=200))
	tanh_label_dict = Dict(1 => "sl_1div80", 2 => "sl_1div90", 3 => "sl_1div100")
	tanh_value_dict = Dict(1 => L÷80, 2 => L÷90, 3 => L÷100)
	# Memory
	l_dict = length(tanh_value_dict)
	data_tanh = zeros(l_dict, 50000, 1024)
	# Time here
	sim_time_tanh = 200:200:10000000
	# Loop over γ
	for i in 1:l_dict
		ggs = zeros(1024)
		gamma_curves_tanh!(ggs, sl=tanh_value_dict[i])
		# Check if there is already a file created
		sim_name_tanh = "gamma_tanh_width_$(tanh_label_dict[i])_tmax_$(sys_tanh.param.Tmax).jld2"
		save_file = string(data_path, sim_name_tanh)
		# If so just read it from disc
		if isfile(save_file)
			println("There is already a file at location $(save_file)\nIf you still want to run the experiment change `sim_name_tanh`, can take several minutes.")
			
		# If not, compute the evolution of the droplet coalescence
		else
	 		data_tanh[i, :, :] = Swalbe.run_gamma(sys_tanh, ggs, r₁=rad, r₂=rad, dump=sys_tanh.param.tdump)
			# Think about storing them on the disc
			df_fluid = Dict()
			# Loop through the time once more
			for t in 1:size(data_tanh)[2]
		        df_fluid["h_$(sim_time_tanh[t])"] = data_tanh[i,t,:]
		    end
			# Save it
			println("Writing file $(tanh_label_dict[i]) to disk")
			save(save_file, df_fluid)
		end
		# Print that you are done
		println("Done with iteration $(tanh_label_dict[i])")
	end
end

# ╔═╡ 0247b632-557c-4151-b89e-c1ebac791205
begin
	plot(fl[1,:], w=3, label="t0")
	plot!(fl[100,:], w=3, label="t100")
	plot!(fl[1000,:], w=3, label="t1000")
	plot!(fl[10000,:], w=3, label="t10000")
end

# ╔═╡ Cell order:
# ╟─c18ba690-9f8e-11ec-1a41-7330ad3642ec
# ╠═a539f57c-7d60-4851-ba5e-60dad313fab7
# ╠═daa88477-ee4d-4563-a0ee-d228a672fd34
# ╠═26520cda-a0ab-403c-9bbe-d39c93d85f7c
# ╟─499e69d3-f894-4c41-b320-021213184c3c
# ╟─cdc87426-d742-4bf0-9699-c793e5753961
# ╟─677a3ed5-046e-49de-a2bf-7f65fe904be4
# ╟─e073abd0-dcd2-4c3d-a043-238f81952446
# ╟─0f007740-d067-4b89-9aeb-2b63f7c77d09
# ╟─f6ebd1b2-fbd3-4622-ab9d-dc98f50bb31f
# ╠═377a8f41-fc08-4941-8b08-530da4f0a987
# ╠═604ce508-e7a6-422f-8a55-e1658143bbac
# ╠═d7e92ccc-a1c5-4547-a6f6-5311514962a1
# ╠═0247b632-557c-4151-b89e-c1ebac791205
