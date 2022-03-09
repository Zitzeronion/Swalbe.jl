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

# ╔═╡ 26520cda-a0ab-403c-9bbe-d39c93d85f7c
include("helpers.jl")

# ╔═╡ c18ba690-9f8e-11ec-1a41-7330ad3642ec
md"# Droplet coalescence - Data generation

This notebook creates all data that we use to understand the coalescence of two sessile droplets on a surface tension gradient.
"

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

# ╔═╡ Cell order:
# ╠═c18ba690-9f8e-11ec-1a41-7330ad3642ec
# ╠═a539f57c-7d60-4851-ba5e-60dad313fab7
# ╠═daa88477-ee4d-4563-a0ee-d228a672fd34
# ╠═26520cda-a0ab-403c-9bbe-d39c93d85f7c
# ╠═499e69d3-f894-4c41-b320-021213184c3c
# ╠═cdc87426-d742-4bf0-9699-c793e5753961
# ╠═677a3ed5-046e-49de-a2bf-7f65fe904be4
# ╠═e073abd0-dcd2-4c3d-a043-238f81952446
