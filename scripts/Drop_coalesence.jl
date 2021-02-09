using DrWatson
@quickactivate :Swalbe
using Plots, CUDA, DataFrames, BSON
CUDA.device!(1)

# Fluid dynamics we need for the experiment
"""
    measure_substratewave

Simulates a changing wettability and the response from the film.
"""
function measure_dropcoalesence(
    sys::Swalbe.SysConst, 
    device::String; 
    R₀ = 1.0, 
    dist = 1,
    dump = 1000,  
    θₛ=ones(sys.Lx, sys.Ly),
    fluid=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    Svelx=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    Svely=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    verbos=true, 
    T=Float64
)
    println("Simulating two coalesing droplets")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    if device == "CPU"
        for i in 1:sys.Lx, j in 1:sys.Ly
            height[i,j] = h₀ + ϵ * sin(2π*wave_x*(i-1)/sys.Lx) * sin(2π*wave_y*(j-1)/sys.Ly)
        end
    elseif device == "GPU"
        h = zeros(size(height))
        for i in 1:sys.Lx, j in 1:sys.Ly
            h[i,j] = h₀ + ϵ * sin(2π*wave_x*(i-1)/sys.Lx) * sin(2π*wave_y*(j-1)/sys.Ly)
        end
        # theta = CUDA.zeros(Float64, sys.Lx, sys.Ly)
        height = CUDA.adapt(CuArray, h)
    end
    Swalbe.equilibrium!(fout, height, velx, vely, vsq)
    ftemp .= fout
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        # Forces are the pressure gradient and the slippage due to substrate liquid boundary conditions
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        # New equilibrium
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        # New moments
        Swalbe.moments!(height, velx, vely, fout)
        # Measurements, in this case only snapshots of simulational arrays
        Swalbe.snapshot!(fluid, height, t, dumping = dump)
        Swalbe.snapshot!(Svelx, velx, t, dumping = dump)
        Swalbe.snapshot!(Svelx, vely, t, dumping = dump)
        Swalbe.snapshot!(theta, θₛ, t, dumping = dump)
        move_substrate!(slipx, θₛ, t, sub_speed, direction=dire)
    end
    return fluid, Svelx, Svely, theta
    CUDA.reclaim()
end

"""
    twodroplets(height, R_drop1, R_drop2, θ, center_drop1, center_drop2)

Initial condition to have two droplets next to each other.
"""
function twodroplets!(height, R_drop1, R_drop2, θ_drop1, θ_drop2, center_drop1, center_drop2)
    # A lot to do here!
    lx, ly = size(height)
    h1 = zeros(lx, ly)
    h2 = zeros(lx, ly)
    for i in 1:lx, j in 1:ly
        circ1 = sqrt((i-center_drop1[1])^2 + (j-center_drop1[2])^2)
        circ2 = sqrt((i-center_drop2[1])^2 + (j-center_drop2[2])^2)
        if circ1 <= R_drop1
            h1[i,j] = (cos(asin(circ1/R_drop1)) - cospi(θ_drop1)) * R_drop1 
        elseif circ2 <= R_drop2
            h1[i,j] += (cos(asin(circ2/R_drop2)) - cospi(θ_drop2)) * R_drop2 
        else
            h1[i,j] = 0.05
        end
    end

    return h1, h2
end

function move_substrate!(θ, input, t, tmove; direction="diagonal")
    if (t % tmove == 0) & (t > 0)
        if direction == "diagonal"
            circshift!(θ, input, (1,1))
        elseif direction == "x"
            circshift!(θ, input, (1,0))
        elseif direction == "y"
            circshift!(θ, input, (0,1))
        end
        input .= θ
    end

    return nothing
end


println("Moving Wettability and possible resonaces")
# Different substrate patches
for direction in ["x"] #  "diagonal"
    # Different initial volumes
    for speed in [0] # 10 100 1000 10000
        for waves in [1 2 3]
            println("Simulating moving substrate wettability with moving direction $(direction) and speed $(speed)")
            sys = Swalbe.SysConst(Lx=512, Ly=512, γ=0.01, δ=1.0, n=3, m=2, hmin=0.07, Tmax=1000000, tdump=1000)
            df_fluid = DataFrame()
            df_velx = DataFrame()
            df_vely = DataFrame()
            df_sub = DataFrame()
            θₚ = ones(sys.Lx,sys.Ly)
            # Substrate patterning
            pattern = "sine"
            if pattern == "sine" 
                for i in 1:sys.Lx, j in 1:sys.Ly
                    θₚ[i,j] = 1/9 + 1/18 * sin(2π*waves*(i-1)/sys.Lx) * sin(2π*waves*(j-1)/sys.Ly)
                end
            elseif pattern == "linear"
                for i in 1:sys.Lx, j in 1:sys.Ly
                    θₚ[i,j] = 1/6 - (i/sys.Lx)/(1/9) 
                end
            end
            # Make a cuarray with the substrate pattern
            θ_in = CUDA.adapt(CuArray, θₚ)
            # Actual simulation
            fluid, velx, vely, substrate = measure_substratewave(sys, "GPU", "no", sub_speed=speed, θₛ=θ_in, dire=direction, dump=sys.tdump)
            println("Writing measurements to dataframes")
            # Filling the dataframes
            for t in 1:sys.Tmax÷sys.tdump
                df_fluid[!, Symbol("h_$(t*sys.tdump)")] = fluid[t,:]
                df_velx[!, Symbol("vx_$(t*sys.tdump)")] = velx[t,:]
                df_vely[!, Symbol("vy_$(t*sys.tdump)")] = vely[t,:]
                df_sub[!, Symbol("theta_$(t*sys.tdump)")] = substrate[t,:]
            end
            println("Saving dataframe subdirection $direction subvel $speed and sines $waves to disk")
            bson("data/Moving_wettability/height_direc_$(direction)_sp_$(speed)_sines_$(waves).bson", Dict(:f => df_fluid))
            bson("data/Moving_wettability/velx_direc_$(direction)_sp_$(speed)_sines_$(waves).bson", Dict(:v => df_velx))
            bson("data/Moving_wettability/vely_direc_$(direction)_sp_$(speed)_sines_$(waves).bson", Dict(:u => df_vely))
            bson("data/Moving_wettability/theta_direc_$(direction)_sp_$(speed)_sines_$(waves).bson", Dict(:t => df_sub))
        
            CUDA.reclaim()
        end
    end
end


println("Script done, let's have a look at the data :)")