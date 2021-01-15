using DrWatson
@quickactivate :Swalbe
using Plots, CUDA, DataFrames, JDF

# Fluid dynamics we need for the experiment
"""
    measure_substratewave

Simulates a changing wettability and the response from the film.
"""
function measure_substratewave(
    sys::Swalbe.SysConst, 
    device::String; 
    h₀ = 1.0, 
    ϵ = 0.1,
    wave_x=1,
    wave_y=1, 
    sub_speed=100,
    dump = 1000,  
    θₛ=ones(sys.Lx, sys.Ly), 
    dire = "x",
    verbos=true, 
    T=Float64
)
    println("Simulating a droplet on a patterned substrate")
    df_fluid = zeros(sys.Tmax÷dump, sys.Lx*sys.Ly) 
    df_v = zeros(sys.Tmax÷dump, sys.Lx*sys.Ly) 
    df_u = zeros(sys.Tmax÷dump, sys.Lx*sys.Ly) 
    df_θ = zeros(sys.Tmax÷sub_speed, sys.Lx*sys.Ly)
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
        theta = CUDA.zeros(sys.Lx, sys.Ly)
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
        Swalbe.snapshot!(df_fluid, height, t, dumping = dump)
        Swalbe.snapshot!(df_v, velx, t, dumping = dump)
        Swalbe.snapshot!(df_u, vely, t, dumping = dump)
        Swalbe.snapshot!(df_θ, θₛ, t, dumping = sub_speed)
        move_substrate!(slipx, θₛ, t, sub_speed, direction=dire)
    end
    return df_fluid, df_v, df_u, df_θ
    CUDA.reclaim()
end

function move_substrate!(θ, input, t, tmove; direction="diagonal")
    if t % tmove == 0 & t > 0
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
for direction in ["x" "diagonal"]
    # Different initial volumes
    for speed in [0 10 100 1000 10000]
        for waves in [1 2 3]
            println("Simulating moving substrate wettability with moving direction $(direction) and speed $(speed)")
            sys = Swalbe.SysConst(Lx=512, Ly=512, γ=0.01, δ=1.0, n=3, m=2, hmin=0.07, Tmax=1000000)
            df_fluid = DataFrame()
            df_velx = DataFrame()
            df_vely = DataFrame()
            df_sub = DataFrame()
            θₚ = ones(sys.Lx,sys.Ly)
            # Substrate patterning
            for i in 1:sys.Lx, j in 1:sys.Ly
                θₚ[i,j] = 1/9 + 1/18 * sin(2π*waves*(i-1)/sys.Lx) * sin(2π*waves*(j-1)/sys.Ly)
            end
            # Make a cuarray with the substrate pattern
            θ_in = CUDA.adapt(CuArray, θₚ)
            # Actual simulation
            fluid, velx, vely, substrate = measure_substratewave(sys, "GPU", sub_speed=speed, θₛ=θ_in, dire=direction, dump=1000)
            println("Writing measurements to dataframes")
            # Filling the dataframes
            for t in 1:sys.Tmax÷sys.tdump
                df_fluid[!, Symbol("h_$(t*sys.tdump)")] = fluid[t,:]
                df_velx[!, Symbol("vx_$(t*sys.tdump)")] = velx[t,:]
                df_vely[!, Symbol("vy_$(t*sys.tdump)")] = vely[t,:]
                df_sub[!, Symbol("theta_$(t*sys.tdump)")] = substrate[t,:]
            end
            println("Saving dataframe subdirection $direction subvel $speed and sines $waves to disk")
            file1 = JDF.save("data/Moving_wettability/height_direc_$(direction)_sp_$(speed)_sines_$(waves).jdf", df_fluid)
            file2 = JDF.save("data/Moving_wettability/velx_direc_$(direction)_sp_$(speed)_sines_$(waves).jdf", df_velx)
            file3 = JDF.save("data/Moving_wettability/vely_direc_$(direction)_sp_$(speed)_sines_$(waves).jdf", df_vely)
            file4 = JDF.save("data/Moving_wettability/theta_direc_$(direction)_sp_$(speed)_sines_$(waves).jdf", df_sub)

            CUDA.reclaim()
        end
    end
end
println("Script done, let's have a look at the data :)")