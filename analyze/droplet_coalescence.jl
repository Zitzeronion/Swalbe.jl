using Plots, DataFrames, JLD2, DataFramesMeta, CSV
using Pkg
Pkg.activate(Base.current_project())
# instantiate, i.e. make sure that all packages are downloaded
Pkg.instantiate()
using Swalbe

#-----------------------------------------------------------#
# 					   System	 							#			
#-----------------------------------------------------------#
L = 1024
γ₀ = 0.0001	
Δγ = 0.00002
TM = 10000000
TLow = 1000
t1000 = 1000:1000:TM
t10 = 10:10:TLow
# Data path
which = "sign"
data_path = "data\\Drop_coalescence_$(which)\\"

#-----------------------------------------------------------#
# 			   Surface tension gardients   					#	
#-----------------------------------------------------------#
function const_gamma(;L=L, γ=γ₀)
    return fill(γ, L)
end
function step_gamma(;L=L, γ=γ₀, perc=20, periodic=false, boundary_l=50)
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
function tanh_gamma(;L=L, γ=γ₀, perc=20, sl=1, periodic=false, boundary_l=50)
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

#-----------------------------------------------------------#
# 					   Run function							#			
#-----------------------------------------------------------#
"""
    run_gamma_periodic(sys, gamma; r₁=115, r₂=115, θ₀=1/9, verbos=true, dump = 100, fluid=zeros(sys.param.Tmax÷dump, sys.L))

Simulation of coalescing droplets on a bounded domain
"""
function run_gamma_periodic(
    sys::Swalbe.SysConst_1D,
    gamma::Vector;
    r₁=115,
    r₂=115, 
    θ₀=1/9,
	drop_cent=(sys.L/3, 2*sys.L/3),  
    verbos=true, 
    dump = 100, 
    fluid=zeros(sys.param.Tmax÷dump, sys.L)
)
    println("Simulating droplet coalecense with surface tension gardient on a periodic domain")
    state = Swalbe.Sys(sys, kind="gamma")
    state.basestate.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state, sys)
    state.γ .= gamma
	Swalbe.∇γ!(state)
	state.∇γ[L-4:L] .= 0
	state.∇γ[1:4] .= 0
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.param.Tmax
        if verbos
            if t % sys.param.tdump == 0
                mass = 0.0
                mass = round(sum(state.basestate.height), digits=3)
				if 0.0 / mass ≠ 0.0
					println("Something went wrong, the mass is likely not a number")
					break
				end
                println("Time step $t bridge height is $(round(minimum(state.basestate.height[sys.L÷2-50:sys.L÷2+50]), digits=3)) and total mass $(mass)")
            end
        end
        Swalbe.filmpressure!(state, sys, γ=gamma)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.basestate.F .= -state.basestate.h∇p .- state.basestate.slip .+ state.∇γ
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(fluid, state.basestate.height, t, dumping = dump)
    end
	
    return fluid
    
end

#-----------------------------------------------------------#
# 					   Plot function						#
#-----------------------------------------------------------#
"""
	do_gif(data_paths; fps=10, end_t=10000)

Generates a gif of droplet coalescence
"""
function do_gif(data_specs; fps=50, end_t=10000, delta_t=1000, folder="Drop_coalescence_sign", L=1024, file_name="default", kind="wall")
	#Check if there is one data set or more
	number_sets = size(data_specs)[1]
	# Have the temporal resolution of the problem
	real_end_T = 1000000*data_specs[1][2]
	time_steps = delta_t:delta_t:real_end_T
	# load data into a DataFrame
	if length(data_specs[1]) == 3
		df = load("data\\$(folder)\\gamma_$(data_specs[1][1])_$(kind)_tmax_$(Int(1000000*data_specs[1][2]))_slip_$(data_specs[1][3])_L_$(L).jld2") |> DataFrame
		# Build the animation object
		drops = @animate for i in 1:fps:(real_end_T÷delta_t)
			time_symbold = Symbol("h_$(time_steps[i])")
			plot(df[!, time_symbold], l=(4, :solid), label="$(data_specs[1][1])")
			if number_sets > 1
				for j in 2:number_sets
					df2 = load("data\\$(folder)\\gamma_$(data_specs[j][1])_$(kind)_tmax_$(Int(1000000*data_specs[1][2]))_slip_$(data_specs[j][3])_L_$(L).jld2") |> DataFrame
					plot!(df2[!, time_symbold], l=(4, :auto), label="$(data_specs[j][1])")
				end
			end
		end
	else
		df = load("data\\$(folder)\\gamma_$(data_specs[1][1])_$(kind)_tmax_$(Int(1000000*data_specs[1][2]))_slip_$(data_specs[1][3])_L_$(L)_hm_$(data_specs[1][4]).jld2") |> DataFrame
		# Build the animation object
		drops = @animate for i in 1:fps:(real_end_T÷delta_t)
			time_symbold = Symbol("h_$(time_steps[i])")
			plot(df[!, time_symbold], l=(4, :solid), label="$(data_specs[1][1])")
			if number_sets > 1
				for j in 2:number_sets
					df2 = load("data\\$(folder)\\gamma_$(data_specs[j][1])_$(kind)_tmax_$(Int(1000000*data_specs[1][2]))_slip_$(data_specs[j][3])_L_$(L)_hm_$(data_specs[j][4]).jld2") |> DataFrame
					plot!(df2[!, time_symbold], l=(4, :auto), label="$(data_specs[j][1])")
				end
			end
		end
	end
	gif(drops, "figures\\coalescence_$(file_name).gif")
end

 #-----------------------------------------------------------#
# 					   Simulations  						#	
# 		  Loop over different configurations				#
#-----------------------------------------------------------#
gamgrads = [const_gamma(),
step_gamma(), 
tanh_gamma(sl=1), 
tanh_gamma(sl=2), 
tanh_gamma(sl=5), 
tanh_gamma(sl=10), 
tanh_gamma(sl=20), 
tanh_gamma(sl=50), 
tanh_gamma(sl=100), 
tanh_gamma(sl=200)]

gamnames = ["const", "step", "tanh1", "tanh2", "tanh5", "tanh10", "tanh20", "tanh50", "tanh100", "tanh200"]

#-----------------------------------------------------------#
# 					Simulation Loop  						#
#-----------------------------------------------------------#
for i in enumerate(gamgrads)
	for j in [15.0]
		sys_loop = Swalbe.SysConstWithBound_1D{Float64}(obs=obst, L=L, 
				param=Swalbe.Taumucs(Tmax=TM, n=9, m=3, δ=j))
		result = run_gamma_walls(sys_loop, i[2], r₁=500, r₂=500, dump=1000)
		data_sim = "gamma_$(gamnames[i[1]])_wall_tmax_$(sys_loop.param.Tmax)_slip_$(Int(sys_loop.param.δ)).jld2"
		save_file = string(data_path, data_sim)
		df_fluid = Dict()
		# Loop through the time once more
		for t in 1:size(result)[1]
		    df_fluid["h_$(t1000[t])"] = result[t,:]
		end
		save(save_file, df_fluid)
		println("Done with simulation $(gamnames[i[1]])")
	end
end
# periodic
for i in enumerate(gamgrads[1:2])
	sll = 11.0
	for j in [0.14] 
		sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=TM, hmin=j, n=9, m=3, δ=sll))
		# sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=TLow, hmin=j, n=9, m=3, δ=sll))
		result = run_gamma_periodic(sys_loop, i[2], r₁=500, r₂=500, dump=1000)
		# result = run_gamma_periodic(sys_loop, i[2], r₁=500, r₂=500, dump=10)
		data_sim = "gamma_$(gamnames[i[1]])_periodic_tmax_$(sys_loop.param.Tmax)_slip_$(Int(sys_loop.param.δ))_L_$(sys_loop.L)_hm_$(Int(round(100*j, digits=2))).jld2"
		save_file = string(data_path, data_sim)
		df_fluid = Dict()
		# Loop through the time once more
		for t in 1:size(result)[1]
		    df_fluid["h_$(t1000[t])"] = result[t,:]
		    # df_fluid["h_$(t10[t])"] = result[t,:]
		end
		save(save_file, df_fluid)
		println("Done with simulation $(gamnames[i[1]])")
	end
end

#-----------------------------------------------------------#
# 				Test loop with vars   						#
#-----------------------------------------------------------#
function do_scan()
	for k in [12.0] 
		for j in [0.12] 
			for i in enumerate(gamgrads[1:2])
				sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=TM, hmin=j, hcrit=0.03, n=9, m=3, δ=k))
				# sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=TLow, hmin=j, n=9, m=3, δ=sll))
				result = run_gamma_periodic_slipcr(sys_loop, i[2], r₁=500, r₂=500, dump=1000)
				# result = run_gamma_periodic(sys_loop, i[2], r₁=500, r₂=500, dump=10)
				data_sim = "gamma_$(gamnames[i[1]])_periodic_tmax_$(sys_loop.param.Tmax)_slip_$(Int(sys_loop.param.δ))_L_$(sys_loop.L)_hm_$(Int(round(100*sys_loop.param.hmin)))_hc_$(Int(round(100*sys_loop.param.hcrit, digits=2))).jld2"
				save_file = string(data_path, data_sim)
				df_fluid = Dict()
				# Loop through the time once more
				for t in 1:size(result)[1]
					df_fluid["h_$(t1000[t])"] = result[t,:]
					# df_fluid["h_$(t10[t])"] = result[t,:]
				end
				save(save_file, df_fluid)
				drops = @animate for i in 1:100:size(result)[1]
					plot(result[i,:])
				end
				gif(drops, "figures\\$(gamnames[i[1]])_slip_$(Int(sys_loop.param.δ))_hm_$(Int(round(100*sys_loop.param.hmin)))_hc_$(Int(round(100*sys_loop.param.hcrit, digits=2))).gif")
				println("Done with simulation $(gamnames[i[1]])")
			end
		end
	end
end
function do_step_scan()
	# TODO: 18.5 last five surface tension gradients starting with tanh10
	gamnames = ["tanh10", "tanh20", "tanh50", "tanh100", "tanh200"] # "const", "step", "tanh1", "tanh2", "tanh5",
	powers = [(9, 3)]
	slips = [12.0]
	hcrits = [0.03]
	hmins = [0.12]
	gammas = [1e-5]
	count = 0
	tHere = 10000:10000:100000000
	for slip in slips
		for k in powers 
			for j in hmins 
				for l in hcrits
					for s in gammas
						# Different surface tension gradients
						gamgrads = [tanh_gamma(sl=10, γ=s), tanh_gamma(sl=20, γ=s), tanh_gamma(sl=50, γ=s), tanh_gamma(sl=100, γ=s), tanh_gamma(sl=200, γ=s)] # const_gamma(γ=s), step_gamma(γ=s), tanh_gamma(sl=1, γ=s), tanh_gamma(sl=2, γ=s), tanh_gamma(sl=5, γ=s), 
						# Loop over gradients
						for gam in enumerate(gamgrads)
							# The system specific constants
							sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=tHere[end], hmin=j, hcrit=l, n=k[1], m=k[2], δ=slip))
							# The actual simulations
							result = run_gamma_periodic(sys_loop, gam[2], r₁=500, r₂=500, dump=tHere[1])
							# Data paths
							data_sim = "gamma_$(Int(round(1000000*s, digits=2)))_$(gamnames[gam[1]])_periodic_tmax_$(sys_loop.param.Tmax)_slip_$(Int(sys_loop.param.δ))_L_$(sys_loop.L)_hm_$(Int(round(100*sys_loop.param.hmin)))_hc_$(Int(round(100*sys_loop.param.hcrit, digits=2)))_gamma_$(Int(round(1000000*s, digits=2))).jld2"
							save_file = string(data_path, data_sim)
							df_fluid = Dict()
							# Loop through the time once more to store the data in a dictonary
							for t in 1:size(result)[1]
								df_fluid["h_$(tHere[t])"] = result[t,:]
							end
							# Save it to disc
							save(save_file, df_fluid)
							# Create an animation
							drops = @animate for i in 1:100:size(result)[1]
								plot(result[i,:])
							end
							gif(drops, "figures\\$(gamnames[gam[1]])_slip_$(Int(sys_loop.param.δ))_hm_$(Int(round(100*sys_loop.param.hmin)))_hc_$(Int(round(100*sys_loop.param.hcrit, digits=2)))_$(k[1]+k[2])_gamma_$(Int(round(1000000*s, digits=2))).gif")
							count += 1
							println("Done with simulation $(count) of $(length(powers)*length(hmins)*length(hcrits)*length(gammas)*length(slips)*length(gamgrads))")
						end
					end
				end
			end
		end
	end
end
#-----------------------------------------------------------#
# 				Data saved on disk   						#
#-----------------------------------------------------------#
#df = load(data_loc) |> DataFrame

#plot(df.h_100, label="t=1", m=(8,:auto), line = (:auto, 4), xlabel="x", ylabel="h", legendfontsize=14, guidefontsize = 18, tickfontsize = 12, legend=true)
#plot!(df.h_100000, label="t=5", m=(8,:auto), line = (:auto, 4))
#plot!(df.h_1000000, label="t=10", m=(8,:auto), line = (:auto, 4))
#plot!(df.h_5000000, label="t=100", m=(8,:auto), line = (:auto, 4))


#-----------------------------------------------------------#
# 					Analysis & Statistics					#			
#-----------------------------------------------------------#
"""
	mymean(h)

First momentum of the thickness field interpreted as a distribution.
"""
function mymean(h)
	μ=0
	m=sum(h)
	for i in 1:length(h)
		μ+= (h[i]/m)*i
	end
	return μ
end

"""
	mystd(h)

Standard deviation of the film thickness.
"""
function mystd(h)
	σ=0
	m=sum(h)
	μ=mymean(h)
	for i in 1:length(h)
		σ += (h[i]/m)*((i-μ)^2)
	end
	return sqrt(σ)
end

"""
	myskew(h)

Computes the skewness, given a thickness field.
"""
function myskew(h)
	σ=mystd(h)
	μ=mymean(h)
	m=sum(h)
	s=0
	for i in 1:length(h)
		s += (h[i]/m)*(((i-μ)/σ)^3)
	end
	return s
end

"""
	bridge_height(df; time=200:200:10000000, L=L, r0=171)

Measurement of bridge height, neck position and skewness.
"""
function bridge_height(df::DataFrame; time=200:200:10000000, hmin=0.11, L=1024, r0=171, slip=12, label_="γ")
	df_ = DataFrame()
	# Parameter
	center = L÷2
	# Controll values
	time_list = []
	grad_list = []
	slip_list = []
	hmin_list = []
	# Computed data
	height_list = []
	skew_list = []
	pos_min_list = []
	# Loop through data
	# for i in values(label_dict)
	# tmp = @subset(df, :sw .== i) 
	for t in time
		neck_region = df[!, Symbol("h_$(t)")][center-r0:center+r0]
		push!(time_list, t)
		push!(grad_list, label_)
		push!(height_list, minimum(neck_region))
		push!(pos_min_list, argmin(neck_region))
		push!(skew_list, myskew(neck_region))
		push!(slip_list, slip)
		push!(hmin_list, hmin)
		# end
	end
	df_[!, "width"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_height"] = height_list
	df_[!, "neck_min"] = pos_min_list
	df_[!, "skewness"] = skew_list
	df_[!, "slip"] = slip_list
	df_[!, "Pi_hmin"] = hmin_list

	return df_
end

"""
	bridge_height(; ; folder="Drop_coalescence_sign", time=200:200:10000000, kind="periodic", slip=12, L=1024, r0=171, label_="step")

Measurement of bridge height, neck position and skewness.
"""
function bridge_height(; folder="Drop_coalescence_sign", T=TM, T2=TLow, hmin=0.11, kind="periodic", slip=12, L=1024, r0=171, label_="step")
	df_ = DataFrame()
	# Parameter
	center = L÷2
	time_dict = Dict(T2 => 10:10:990, T => 1000:1000:TM)
	# Controll values
	time_list = []
	grad_list = []
	hm_list = []
	slip_list = []
	# Computed data
	min_list = []
	max_list = []
	skew_list = []
	pos_min_list = []
	pos_max_list = []
	# Load data
	for i in [T2, T]
		df = load("data\\$(folder)\\gamma_$(label_)_$(kind)_tmax_$(i)_slip_$(slip)_L_$(L)_hm_$(Int(100*hmin)).jld2") |> DataFrame
		println("Reading file:\ndata\\$(folder)\\gamma_$(label_)_$(kind)_tmax_$(i)_slip_$(slip)_L_$(L)_hm_$(Int(100*hmin)).jld2 ")
		# Extract data from the simulation
		for t in time_dict[i]
			neck_region = df[!, Symbol("h_$(t)")][center-r0:center+r0]
			push!(time_list, t)
			push!(grad_list, label_)
			push!(min_list, minimum(neck_region))
			push!(max_list, maximum(neck_region))
			push!(pos_min_list, argmin(neck_region))
			push!(pos_max_list, argmax(neck_region))
			push!(skew_list, myskew(neck_region))
			push!(slip_list, slip)
			push!(hm_list, hmin)
		end
	end
	# Fill the dataframe with the results of the analysis
	df_[!, "gamma"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_min"] = min_list
	df_[!, "bridge_max"] = max_list
	df_[!, "neck_pos"] = pos_min_list
	df_[!, "neck_max_pos"] = pos_max_list
	df_[!, "skewness"] = skew_list
	df_[!, "slip"] = slip_list
	df_[!, "hmin"] = hm_list
	# Return the denser data
	return df_
end

function data_loop(labels; slip=12, extras="hm_11")
	df = DataFrame()
	for i in labels
		tmp = bridge_height(; label_=i, slip=slip)
		df = vcat(df, tmp)
		println("Done with data set: $(i)")
	end
	CSV.write("data\\coalescence_analysed\\bridge_data_$(T)_slip_$(slip)_$(extras).csv", df)
	return df
end

df_drop = CSV.read("data\\coalescence_analysed\\bridge_data.csv", DataFrame)

time_set = 100:100:5000000
data_const = bridge_height(df, time=100:100:5000000, label_="γ₀")
plot(time_set, 
     data_const.bridge_height, 
     label="γ₀",
     axis=:log, 
     m=(8,:circle), 
     line = (:auto, 4), 
     xlabel="t", 
     ylabel="h", 
     legendfontsize=14, 
     guidefontsize = 18, 
     tickfontsize = 12, 
     legend=:bottomright)
plot!(time_set, 0.001 .* time_set.^(2/3), label="t^(2/3)")