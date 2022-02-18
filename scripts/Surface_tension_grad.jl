using Plots, DataFrames, FileIO
using Swalbe
# Constants
L = 1024
x = collect(1:L)
γ = zeros(4,L)
ε = 0.02
γ₀ = 0.0001
γ_bar = (γ₀ + (γ₀ - ε))/2
Δγ = ε
sl = L÷10

# Surface tension functions
function gamma_curves!(x; x0=γ₀, ϵ=0.02, L=L, sl=sl)
	function smooth(x, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((x .- L÷2) ./ (sl))))
	end
	x[1,:] .= x0
	x[2,:] .= x0 .* (1 .- ϵ .* collect(1:L) ./ L)
	x[3,1:L÷2] .= x0
	x[3,L÷2+1:L] .= x0 - x0 * ϵ
	# x[3,:] .= x0 .* exp.(-collect(1:L)/L)
	x[4,:] .= x0 .* smooth(x[3,:], L, sl) .+ (1 .- smooth(x[3,:], L, sl)) .* x0 .*(1 - ϵ) 
end

# Initial state
rad = 500
h = Swalbe.two_droplets(Swalbe.SysConst_1D(L=1024), r₁=rad, r₂=rad)
p0 = plot(collect(1:1024), h, 
		  w=3, 
		  aspect_ratio=7, 
		  label="Initial conf.", 
		  xlabel="x [Δx]", 
		  ylabel="h(x)",
		  legendfontsize = 14,			# legend font size
          tickfont = (14),	# tick font and size
          guidefont = (15))
ylims!(0,40)

# Run function to perform the experiments
function run_(
    sys::Swalbe.SysConst_1D,
    gamma::Vector;
    r₁=115,
    r₂=115, 
    θ₀=1/9,  
    verbos=true, 
    dump = 100, 
    fluid=zeros(sys.Tmax÷dump, sys.L)
)
    println("Simulating droplet coalecense with surface tension gardient")
    state = Swalbe.Sys(sys, kind="gamma")
    drop_cent = (sys.L/3, 2*sys.L/3)
    state.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state)
    state.γ .= gamma
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
        Swalbe.∇γ!(state)
        state.F .= -state.h∇p .- state.slip .- state.γ
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(fluid, state.height, t, dumping = dump)
    end

    return fluid
    
end

# Define a SysConst and try if the simulation runs
sys = Swalbe.SysConst_1D(L=1024, Tmax=4000000, δ=50.0)
gamma_curves!(γ, x0=1e-4)
l = run_(sys, γ[2,:], r₁=rad, r₂=rad)

# Loop through the different surface tensions and disjoining pressure terms
data = zeros(40000, L, 8)
for i in 1:8
    k = 0
    if i < 5
        sys = Swalbe.SysConst_1D(L=1024, Tmax=4000000, δ=50.0)
        k = i
    else
        sys = Swalbe.SysConst_1D(L=1024, n=3, m=2, Tmax=4000000, δ=50.0)
        k = i - 4
    end 
    data[:,:,k] = run_(sys, γ[k,:], r₁=rad, r₂=rad)
    println("Done with iteration $i")
end

p1 = plot(data[100,:,1], label="t100"); plot!(data[10000, :, 1], label="1e6"); plot!(data[40000, :, 1], label="t4e6")
p2 = plot(data[100,:,2], label="t100"), plot!(data[10000, :, 2], label="1e6"); plot!(data[40000, :, 2], label="t4e6")
p3 = plot(data[100,:,3], label="t100"); plot!(data[10000, :, 3], label="1e6"); plot!(data[40000, :, 3], label="t4e6")
p4 = plot(data[100,:,4], label="t100"); plot!(data[10000, :, 4], label="1e6"); plot!(data[40000, :, 4], label="t4e6")

plot(p1, p2, p3, p4, layout = (2, 2), legend = false)