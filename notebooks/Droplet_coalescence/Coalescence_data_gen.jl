### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ c18ba690-9f8e-11ec-1a41-7330ad3642ec
md"# Droplet coalescence - Data generation

This notebook creates all data that we use to understand the coalescence of two sessile droplets on a surface tension gradient.
"

# ╔═╡ c4d277dd-6c3d-4f99-a084-57d213e9e689
begin
	data_path = "..\\..\\data\\Drop_coalescence\\"
	L = 1024
	x = collect(1:L)
	γ = zeros(4,L)
	ε = 0.2
	γ₀ = 0.0001
	γ_bar = (γ₀ + (γ₀ - ε))/2
	Δγ = ε
	sl = L÷10
	gamma_labels = Dict(1 => "default", 2 => "linear", 3 => "step", 4 => "tanh")
end

# ╔═╡ 18b7d0c7-dbe0-4afe-9afb-bef6e6437d43
"""
	gamma_curves!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)

Different surface tension fields, constant, linear, step and tanh.
"""
function gamma_curves!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	# Smoothing
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	# Constant surface tension
	out[1,:] .= x0
	# Linear depency
	out[2,:] .= x0 .* (1 .- ϵ .* l ./ L)
	# Step function
	out[3,1:L÷2] .= x0
	out[3,L÷2+1:L] .= x0 - x0 * ϵ
	# Tangent hyperbolicus smoothing
	out[4,:] .= x0 .* smooth(x, L, sl) .+ (1 .- smooth(x, L, sl)) .* x0 .*(1 - ϵ) 
	return nothing
end

# ╔═╡ 9d2c1769-4c93-4926-b5bb-62dc577fce6e
"""
	gamma_curves_tanh!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)

Different surface tension fields using the tangent hyperbolicus and a varying smoothing width.
"""
function gamma_curves_tanh!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	
	out[:] .= x0 .* smooth(x, L, sl) .+ (1 .- smooth(x, L, sl)) .* x0 .*(1 - ϵ) 
	return nothing
end

# ╔═╡ ecbce468-bdc8-40df-b8dd-d2b1bb7ffd01
"""
	run_(sys::Swalbe.SysConst_1D,
    	 gamma::Vector;
    	 r₁=115,
    	 r₂=115, 
    	 θ₀=1/9,  
    	 verbos=true, 
    	 dump = 100, 
    	 fluid=zeros(sys.param.Tmax÷dump, sys.L))

Lattice Boltzmann simulation of coalescing droplets.
"""
function run_(
    sys::Swalbe.SysConst_1D,
    gamma::Vector;
    r₁=115,
    r₂=115, 
    θ₀=1/9,  
    verbos=true, 
    dump = 100, 
    fluid=zeros(sys.param.Tmax÷dump, sys.L)
)
    println("Simulating droplet coalecense with surface tension gardient")
    state = Swalbe.Sys(sys, kind="gamma")
    drop_cent = (sys.L/3, 2*sys.L/3)
    state.basestate.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state, sys)
    state.γ .= gamma
	Swalbe.∇γ!(state)
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.basestate.height)
            if verbos
                println("Time step $t bridge height is $(round(minimum(state.basestate.height[sys.L÷2-20:sys.L÷2+20]), digits=3))")
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
	# println("Max γ: $(maximum(state.γ)),\nMin γ: $(minimum(state.γ))")
    return fluid
    
end

# ╔═╡ 499e69d3-f894-4c41-b320-021213184c3c
md"## First run

Up to 5×10⁶ time step with the four surface tension curves."

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
	 		data[i, :, :] = run_(sys, γ[i, :], r₁=rad, r₂=rad)
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

Up to 10⁷ time step with different smoothing width in the tangent hyperbolicus."

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
			runyesno = readline()
		# If not, compute the evolution of the droplet coalescence
		else
	 		data_tanh[i, :, :] = run_(sys_tanh, ggs, r₁=rad, r₂=rad, dump=sys_tanh.param.tdump)
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

# ╔═╡ Cell order:
# ╠═c18ba690-9f8e-11ec-1a41-7330ad3642ec
# ╠═a539f57c-7d60-4851-ba5e-60dad313fab7
# ╠═daa88477-ee4d-4563-a0ee-d228a672fd34
# ╠═c4d277dd-6c3d-4f99-a084-57d213e9e689
# ╠═18b7d0c7-dbe0-4afe-9afb-bef6e6437d43
# ╠═9d2c1769-4c93-4926-b5bb-62dc577fce6e
# ╠═ecbce468-bdc8-40df-b8dd-d2b1bb7ffd01
# ╠═499e69d3-f894-4c41-b320-021213184c3c
# ╠═cdc87426-d742-4bf0-9699-c793e5753961
# ╠═677a3ed5-046e-49de-a2bf-7f65fe904be4
# ╠═e073abd0-dcd2-4c3d-a043-238f81952446
