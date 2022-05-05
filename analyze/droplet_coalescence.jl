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
# Boundaries
obst = zeros(L)
obst[1]=1
obst[2]=1
obst[end]=1
obst[end-1]=1
# Data path
which = "sign"
data_path = "data\\Drop_coalescence_$(which)\\"

#-----------------------------------------------------------#
# 			   Surface tension gardients   					#	
#-----------------------------------------------------------#
function const_gamma(;L=L, γ=γ₀)
    return fill(γ, L)
end
function step_gamma(;L=L, γ=γ₀, Δ=Δγ, periodic=false, boundary_l=50)
    x = ones(L)
    for i in 1:L
        if i < L÷2
            x[i] = γ
        else
            x[i] = γ - Δ
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
function tanh_gamma(;L=L, γ=γ₀, Δ=Δγ, sl=1, periodic=false, boundary_l=50)
    l = collect(1.0:L)
    function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	x = ones(L)
	x .= γ .* smooth(l, L, sl) .+ (1 .- smooth(l, L, sl)) .* (γ  - Δ) 
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
    run_gamma_walls(sys, gamma; r₁=115, r₂=115, θ₀=1/9, verbos=true, dump = 100, fluid=zeros(sys.param.Tmax÷dump, sys.L))

Simulation of coalescing droplets on a bounded domain
"""
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
                mass = round(sum(state.basestate.height), digits=3)
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

"""
    run_gamma_nograd(sys, gamma; r₁=115, r₂=115, θ₀=1/9, verbos=true, dump = 100, fluid=zeros(sys.param.Tmax÷dump, sys.L))

Simulation of coalescing droplets on a bounded domain
"""
function run_gamma_nograd(
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
    println("Simulating droplet coalecense with surface tension gardient on a domain with walls")
    state = Swalbe.Sys(sys, kind="gamma")
    state.basestate.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state, sys)
    state.γ .= gamma
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.param.Tmax
        if verbos
            if t % sys.param.tdump == 0
                mass = 0.0
                mass = round(sum(state.basestate.height), digits=3)
                println("Time step $t bridge height is $(round(minimum(state.basestate.height[sys.L÷2-50:sys.L÷2+50]), digits=3)) and total mass $(mass)")
            end
        end
        Swalbe.filmpressure!(state, sys, γ=gamma)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.basestate.F .= -state.basestate.h∇p .- state.basestate.slip #.+ state.∇γ
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
	four_plots(f; t1=1, t2=100, t3=1000, t4=10000)

Plots four iterations of the coalescence simulation
"""
function four_plots(f; t1=1, t2=100, t3=1000, t4=10000)
    plot(f[t1, :], xlabel="x", ylabel="h", l=(4, :auto), label="t10^2")
    plot!(f[t2, :], l=(4, :auto), label="t10^3")
    plot!(f[t3, :], l=(4, :auto), label="t10^4")
    plot!(f[t4, :], l=(4, :auto), label="t10^5")
end

"""
	do_gif(data_paths; fps=10, end_t=10000)

Generates a gif of droplet coalescence
"""
function do_gif(data_specs; fps=50, end_t=10000, delta_t=1000, folder="Drop_coalescence_sign", L=1024, file_name="default", kind="wall")
	#Check if there is one data set or more
	number_sets = size(data_specs)[1]
	# Have the temporal resolution of the problem
	time_steps = delta_t:delta_t:end_t
	# load data into a DataFrame
	df = load("data\\$(folder)\\gamma_$(data_specs[1][1])_$(kind)_tmax_$(data_specs[1][2])_slip_$(data_specs[1][3])_L_$(L).jld2") |> DataFrame
	# Build the animation object
	drops = @animate for i in 1:fps:(end_t÷delta_t)
		time_symbold = Symbol("h_$(time_steps[i])")
		plot(df[!, time_symbold], l=(4, :solid), label="$(data_specs[1][1])")
		if number_sets > 1
			for j in 2:number_sets
				df2 = load("data\\$(folder)\\gamma_$(data_specs[j][1])_$(kind)_tmax_$(data_specs[j][2])_slip_$(data_specs[j][3])_L_$(L).jld2") |> DataFrame
				plot!(df2[!, time_symbold], l=(4, :auto), label="$(data_specs[j][1])")
			end
		end
	end
	gif(drops, "figures\\coalescence_$(file_name).gif")
end


#-----------------------------------------------------------#
# 					   Simulation   						#	
# 			Constant surface tension plain coal.			#
#-----------------------------------------------------------#
sys_bound = Swalbe.SysConstWithBound_1D{Float64}(obs=obst, L=L, 
			param=Swalbe.Taumucs(Tmax=1000000, n=9, m=3, δ=15.0, hcrit=0.05, hmin=0.098))
fl = run_gamma_walls(sys_bound, tanh_gamma(sl=10), r₁=500, r₂=500, dump=1000)
# How does it look like
four_plots(fl)
# How does the bridge look like
bmin = []
for i in 1:10000
    push!(bmin, minimum(fl[i, 250:750]))
end
t100 = 100:100:1000000 
plot(t100,bmin, xlabel="t", axis=:log, ylabel="h_b", l=(4, :auto), label="bridge", legend=:bottomright)
plot!(t100, 0.005 .* t100.^(2/3), l=(4, :auto), label="fit", ylim=(0.1, 20))

#-----------------------------------------------------------#
# 					   Simulation   						#	
# 		Surface tension gradient with tanh(sl=10)			#
#-----------------------------------------------------------#
g2 = tanh_gamma(sl=10)
fl_tanh = run_gamma_walls(sys_bound, g2, r₁=500, r₂=500, dump=1000)
four_plots(fl_tanh)


sys_per = Swalbe.SysConst_1D(L=L, 
			param=Swalbe.Taumucs(Tmax=10000, n=9, m=3, δ=12.0))
#-----------------------------------------------------------#
# 					   Simulations  						#	
# 		  Loop over different configurations				#
#-----------------------------------------------------------#
gamgrads_red = [const_gamma(),
			    step_gamma(), 
			    tanh_gamma(sl=50)]
			

gamnames_red = ["const", "step", "tanh50"]
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
for i in enumerate(gamgrads)
	for j in [6.0]
		# sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=TM, n=9, m=3, δ=j))
		sys_loop = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=TLow, n=9, m=3, δ=j))
		# result = run_gamma_periodic(sys_loop, i[2], r₁=500, r₂=500, dump=1000)
		result = run_gamma_periodic(sys_loop, i[2], r₁=500, r₂=500, dump=10)
		data_sim = "gamma_$(gamnames[i[1]])_periodic_tmax_$(sys_loop.param.Tmax)_slip_$(Int(sys_loop.param.δ))_L_$(sys_loop.L).jld2"
		save_file = string(data_path, data_sim)
		df_fluid = Dict()
		# Loop through the time once more
		for t in 1:size(result)[1]
		    # df_fluid["h_$(t1000[t])"] = result[t,:]
		    df_fluid["h_$(t10[t])"] = result[t,:]
		end
		save(save_file, df_fluid)
		println("Done with simulation $(gamnames[i[1]])")
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
function bridge_height(df::DataFrame; time=200:200:10000000, L=1024, r0=171, label_="γ")
	df_ = DataFrame()
	# Parameter
	center = L÷2
	# Controll values
	time_list = []
	grad_list = []
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
		# end
	end
	df_[!, "width"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_height"] = height_list
	df_[!, "neck_min"] = pos_min_list
	df_[!, "skewness"] = skew_list

	return df_
end

"""
	bridge_height(; ; folder="Drop_coalescence_sign", time=200:200:10000000, kind="periodic", slip=12, L=1024, r0=171, label_="step")

Measurement of bridge height, neck position and skewness.
"""
function bridge_height(; folder="Drop_coalescence_sign", T=TM, time=1000:1000:10000000, kind="periodic", slip=12, L=1024, r0=171, label_="step")
	df_ = DataFrame()
	# Parameter
	center = L÷2
	# Controll values
	time_list = []
	grad_list = []
	# Computed data
	min_list = []
	max_list = []
	skew_list = []
	pos_min_list = []
	pos_max_list = []
	# Load data
	df = load("data\\$(folder)\\gamma_$(label_)_$(kind)_tmax_$(T)_slip_$(slip)_L_$(L).jld2") |> DataFrame
	println("Reading file:\ndata\\$(folder)\\gamma_$(label_)_$(kind)_tmax_$(T)_slip_$(slip)_L_$(L).jld2 ")
	# Extract data from the simulation
	for t in time
		neck_region = df[!, Symbol("h_$(t)")][center-r0:center+r0]
		push!(time_list, t)
		push!(grad_list, label_)
		push!(min_list, minimum(neck_region))
		push!(max_list, maximum(neck_region))
		push!(pos_min_list, argmin(neck_region))
		push!(pos_max_list, argmax(neck_region))
		push!(skew_list, myskew(neck_region))
	end
	# Fill the dataframe with the results of the analysis
	df_[!, "gamma"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_min"] = min_list
	df_[!, "bridge_max"] = max_list
	df_[!, "neck_pos"] = pos_min_list
	df_[!, "neck_max_pos"] = pos_max_list
	df_[!, "skewness"] = skew_list
	# Return the denser data
	return df_
end

function data_loop(labels; time=1000:1000:10000000, T=TM, slip=12)
	df = DataFrame()
	for i in labels
		tmp = bridge_height(; label_=i, time=time, T=T, slip=slip)
		df = vcat(df, tmp)
		println("Done with data set: $(i)")
	end
	CSV.write("data\\coalescence_analysed\\bridge_data_$(T)_slip_$(slip).csv", df)
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