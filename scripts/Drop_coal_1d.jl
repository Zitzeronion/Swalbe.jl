using DataFrames, CairoMakie

# Define initial condition
"""
    twodroplets(L)

Generates a fluid state of two droplets with radii r₁, r₂, contact angles θ₁, θ₂ and centers `center`.
"""
function twodroplets(L; r₁=230, r₂=230, θ₁=1/9, θ₂=1/9, center=(L/3,2L/3))
    dum = zeros(L)
    dum2 = zeros(L)
    height = zeros(L)
    # area = 2π * radius^2 * (1- cospi(θ))
    @inbounds for i in 1:L
        circ = sqrt((i-center[1])^2)
        if circ <= r₁
            dum[i] = (cos(asin(circ/r₁)) - cospi(θ₁)) * r₁
        else    
            dum[i] = 0.05
        end
    end
    
    @inbounds for i in 1:L
        circ2 = sqrt((i-center[2])^2)
        if circ2 <= r₂
            dum2[i] = (cos(asin(circ2/r₂)) - cospi(θ₂)) * r₂
        else    
            dum2[i] = 0.05
        end
    end

    @inbounds for i in 1:L
        if dum[i] < 0
            dum[i] = 0.05
        end
        if dum2[i] < 0
            dum2[i] = 0.05
        end
    end

    height .= dum .+ dum2 .- 0.05
    @inbounds for i in 1:L
        if height[i] < 0
            height[i] = 0.05
        end
    end
    return height
end
time_dump = 100
"""
    run_drop_coal()

Simulation for droplet coalesence in quasi two dimensions.
"""
function run_drop_coal(
    sys::Swalbe.SysConst_1D;
    r₁=115,
    r₂=115, 
    θ₀=1/9,  
    verbos=true, 
    dump = time_dump, 
    fluid=zeros(sys.Tmax÷dump, sys.L),
    T=Float64
)
    println("Simulating droplet coalecense")
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p, thermal = Swalbe.Sys(sys, true, T)
    drop_cent = (sys.L/3, 2*sys.L/3)
    height .= twodroplets(sys.L, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(feq, height, vel)
    ftemp .= feq
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t bridge height is $(round(height[Int(sys.L/2)], digits=3))")
            end
        end
        # Film pressure with rho field
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        Swalbe.thermal!(thermal, height, sys.kbt, sys.μ, sys.δ)
        # Here we a force that is like pull of an inclined plane
        F .= h∇p .+ slip .+ thermal 
        Swalbe.equilibrium!(feq, height, vel)
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        Swalbe.moments!(height, vel, fout)
        # Make a snaphot of the configuration
        Swalbe.snapshot!(fluid, height, t, dumping = dump)
    end

    return fluid# height, vel, rho
    
end
# Now run the actual simulation

sys = Swalbe.SysConst_1D(L=1024, Tmax=2000000, δ=50.0, γ=0.005)
fluid = run_drop_coal(sys, r₁=500, r₂=500)

sys2 = Swalbe.SysConst_1D(L=1024, Tmax=2000000, δ=10.0, γ=0.005)
fluid2 = run_drop_coal(sys2, r₁=500, r₂=500)

sys3 = Swalbe.SysConst_1D(L=1024, Tmax=2000000, δ=200.0, γ=0.0001)
fluid3 = run_drop_coal(sys3, r₁=500, r₂=500)

bridge_h = Float64[]
bridge_h2 = Float64[]
bridge_h3 = Float64[]
for i in 1:size(fluid)[1]
    push!(bridge_h, fluid[i, sys.L÷2])
    push!(bridge_h2, fluid2[i, sys.L÷2])
    push!(bridge_h3, fluid3[i, sys.L÷2])
end

let
    lines(1:100:2000000, bridge_h, figure = (resolution = (700,450),),
        axis = (xscale = log10, yscale = log10, xlabel = L"time [Δt]", ylabel = L"h_B",
        xgridstyle=:dash, ygridstyle=:dash, xminorticksvisible = true,
        xminorticks = IntervalsBetween(9), yminorticksvisible = true,
        yminorticks = IntervalsBetween(9)))
    current_figure()
end

lines!(1:100:2000000, bridge_h2)
lines!(1:100:2000000, bridge_h3)
current_figure()



# Push the results into a dict for writting to CSV
df_fluid = Dict()
df_rho = Dict()
# Fill the dict with the time steps recorded
for t in 1:sys.Tmax÷time_dump
    df_fluid["h_$(t*time_dump)"] = fluid[t,:]
end
# Write the desired CSV files
#CSV.write("drop_coal_data/height_run_drop_coal_first_try.csv", df_fluid)
save("drop_coal_data/height_run_drop_coal_thrid_try.jld2", df_fluid)
# Simulation is done and output is save to files height_...csv and rho_...csv
println("Sim done, results in height_{...}.csv")