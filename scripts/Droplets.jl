using Swalbe, Plots


"""
    run_drop_coal()

Simulation of two droplets in close contact to study their coalesence behaviour. 
"""
function run_drop_coal(
    sys::Swalbe.SysConst_1D;                # System relevant constants
    r₁=500,                                 # First droplets radius
    r₂=500,                                 # Second droplets radius
    θ₀=1/9,                                 # Contact angle 
    centers=(sys.L/3, 2*sys.L/3),           # Drople centers 
    dump = 100,                             # Data dumping interval
    fluid=zeros(sys.Tmax÷dump, sys.L),      # Data
    verbos=true                             # Simulation updates about the bridge height in consol
)
    println("Simulating droplet coalecense")
    state = Swalbe.Sys(sys)
    state.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=centers)
    Swalbe.equilibrium!(state)
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbos
                println("Time step $t bridge height is $(round(state.height[Int(sys.L/2)], digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip 
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
        # Make a snaphot of the configuration
        Swalbe.snapshot!(fluid, state.height, t, dumping = dump)
    end
    return fluid
end

"""
    bridge_height(data)

Compuation of the bridge height between two droplets based on input `data`.
"""
function bridge_height(data; t1=100, t_dur=size(data)[1], pos=size(data)[2]÷2)
    bridge_h = []
    for i in 1:t_dur
        push!(bridge_h, data[i, pos])
    end
    return bridge_h
end

"""
	τν(γ, R, η)

Computes the viscous time.

# Mathematics

``t_c = \\sqrt{\\frac{\\rho R^3}{\\gamma}}``
"""
function τν(γ; R=500, η=1/6)
	time_c = 0.0
	time_c = R*η/γ
	return time_c
end

# Different surface tension values
γs = [0.00001, 0.00005, 0.0001, 0.0005]
# Large array that contains the simulation data
data_merge = zeros(20000, 1024, length(γs))
# Loop different surface tension values
for γ in enumerate(γs)
    sys = Swalbe.SysConst_1D(L=1024, Tmax=2000000, δ=50.0, γ=γ[2])
    data_merge[:,:,γ[1]] = run_drop_coal(sys, r₁=sphere_rad, r₂=sphere_rad)
end
# Compute the bridge heights for the data collect and store it in a single array
bridges = zeros(20000, length(γs))
for i in 1:length(γs)
	bridges[:, i] = bridge_height(data_merge[:,:,i])
end
# Actual LBM time, time steps without (physical) meaning
time_lbm = 100:100:2000000
p1 = plot(time_lbm,                                                 # x-data
		  bridges,                                                  # y-data
		  xlabel="t [Δt]",                                          # axis label
		  ylabel="h₀",                                              # axis label
		  label=["γ=0.00001" "γ=0.00005" "γ=0.0001" "γ=0.0005"],    # labels for the curves
		  xticks=([0:1000000:2000000;], ["0", "1×10⁶", "2×10⁶"]),   # ticks and labels
		  xlims=(0,2100000),                                        # well limits of display
		  ylims=(0, 40),                                            # well limits of display
		  w = 3, 							                        # line width
		  marker = (8, :auto, 0.6),			                        # marker size 
		  legend=:bottomright,                                      # legend position
		  legendfontsize = lf,			                            # legend font size
          tickfont = (14),	                                        # tick font and size
          guidefont = (15),	                                        # label font and size
		  grid=:none                                                # no gird, was asked once by PRE to remove the grid
)
# Now to collapse the curves we non-dimensionalize the time with the visous time scale τν and the height with the sphere radius R₀
time_norm = zeros(length(time_lbm), length(γs))
for i in enumerate(γs)
	time_norm[:, i[1]] .= time_lbm ./ τν(i[2])
end
# And plot it again
p2 = plot(time_norm,                                                # x-data
		  bridges ./ 500.0,                                         # y-data
		  xlabel="t [Δt]",                                          # axis label
		  ylabel="h₀",                                              # axis label
		  label=["γ=0.00001" "γ=0.00005" "γ=0.0001" "γ=0.0005"],    # labels for the curves
		  xticks=([0:1000000:2000000;], ["0", "1×10⁶", "2×10⁶"]),   # ticks and labels
		  xlims=(0,2100000),                                        # well limits of display
		  ylims=(0, 40),                                            # well limits of display
		  w = 3, 							                        # line width
		  marker = (8, :auto, 0.6),			                        # marker size 
		  legend=:bottomright,                                      # legend position
		  legendfontsize = lf,			                            # legend font size
          tickfont = (14),	                                        # tick font and size
          guidefont = (15),	                                        # label font and size
		  grid=:none                                                # no gird, was asked once by PRE to remove the grid
)
# In the inertial regime we assume that the bridge height grows as time to the power of 2/3, we can check it by simply adding this curve
t_scaling = collect(1e-5:1e-5:10)                                       # Finer time scale for the power law
plot!(t_scaling, 0.024 .* t_scaling.^(2/3), c=:black, w=4, l=:dash, label="t^(2/3)")
