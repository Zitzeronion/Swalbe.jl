using DrWatson
@quickactivate :Swalbe
using CUDA, DataFrames, BSON
# CUDA.device!(1)

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
    fluid=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    theta=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    dire = "x",
    verbos=true, 
    T=Float64
)
    println("Simulating a droplet on a patterned substrate")
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
        Swalbe.snapshot!(theta, θₛ, t, dumping = dump)
        move_substrate!(slipx, θₛ, t, sub_speed, direction=dire)
    end
    return fluid, theta
    CUDA.reclaim()
end
# Well not most clean way but one one to freeze the substrate without having the x/0 issue.
function measure_substratewave(
    sys::Swalbe.SysConst, 
    device::String,
    move_sub::String;
    h₀ = 1.0, 
    ϵ = 0.1,
    wave_x=1,
    wave_y=1, 
    sub_speed=100,
    dump=1000,  
    θₛ=ones(sys.Lx, sys.Ly),
    fluid=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    theta=zeros(sys.Tmax÷dump, sys.Lx*sys.Ly),
    dire = "x",
    verbos=true, 
    T=Float64
)
    println("Simulating a time dependent substrate pattern")
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
        Swalbe.snapshot!(theta, θₛ, t, dumping = dump)
        if move_sub == "yes"
            # move_substrate!(slipx, θₛ, t, sub_speed, direction=dire)
        end
    end
    return fluid, theta
    CUDA.reclaim()
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

function pyramidpattern(Lx, Ly; waves=1, mid=1/9, lower=1/18)
    half = Lx÷2
    quat = Lx÷4
    sixt = Lx÷6
    space = zeros(Lx,Ly)
    theta = zeros(Lx,Ly)
    if waves == 1
        for i in 1:Lx÷2
            space[i,i:(half+1-i)] .= (i-1)
            space[half+1-i,i:(half+1-i)] .= (i-1)
            space[i:(half+1-i),i] .= (i-1)
            space[i:(half+1-i),half+1-i] .= (i-1)
        end
        space[half+1:Lx, 1:half] .= - space[1:half,1:half]
        space[1:half, half+1:Lx] .= - space[1:half,1:half]
        space[half+1:Lx, half+1:Lx] .= space[1:half,1:half]
        theta .= mid .+ lower .* space ./ half
    
    elseif waves == 2
        for i in 1:Lx÷4
            space[i,i:(quat+1-i)] .= (i-1)
            space[quat+1-i,i:(quat+1-i)] .= (i-1)
            space[i:(quat+1-i),i] .= (i-1)
            space[i:(quat+1-i),quat+1-i] .= (i-1)
        end
        space[quat+1:half, 1:quat] .= - space[1:quat,1:quat]
        space[1:quat, quat+1:half] .= - space[1:quat,1:quat]
        space[quat+1:half, quat+1:half] .= space[1:quat,1:quat]

        space[1:half, half+1:Lx] .= space[1:half,1:half]
        space[half+1:Lx, half+1:Lx] .= space[1:half,1:half]
        space[half+1:Lx, 1:half] .= space[1:half,1:half]
        theta .= mid .+ lower .* space ./ quat

    elseif waves == 3
        thr = Lx÷3
        for i in 1:Lx÷6
            space[i,i:(sixt+1-i)] .= (i-1)
            space[sixt+1-i,i:(sixt+1-i)] .= (i-1)
            space[i:(sixt+1-i),i] .= (i-1)
            space[i:(sixt+1-i),sixt+1-i] .= (i-1)
        end
        space[sixt+1:thr, 1:sixt] .= - space[1:sixt,1:sixt]
        space[1:sixt, sixt+1:thr] .= - space[1:sixt,1:sixt]
        space[sixt+1:thr, sixt+1:thr] .= space[1:sixt,1:sixt]
        
        space[1:thr, thr+1:2*thr] .= space[1:thr,1:thr]
        space[1:thr, 2*thr+1:3*thr] .= space[1:thr,1:thr]
        space[thr+1:2*thr, 1:thr] .= space[1:thr,1:thr]
        space[thr+1:2*thr, thr+1:2*thr] .= space[1:thr,1:thr]
        space[thr+1:2*thr, 2*thr+1:3*thr] .= space[1:thr,1:thr]
        space[2*thr+1:3*thr, 1:thr] .= space[1:thr,1:thr]
        space[2*thr+1:3*thr, thr+1:2*thr] .= space[1:thr,1:thr]
        space[2*thr+1:3*thr, 2*thr+1:3*thr] .= space[1:thr,1:thr]
        theta .= mid .+ lower .* space ./ sixt
    end

    return theta
end

println("Moving Wettability and possible resonaces")
v_lam1_dia = [0, 4900, 490, 49]
v_lam2_dia = [0, 9802, 980, 98]
v_lam3_dia = [0, 14702, 1470, 147] 
for direction in ["diagonal"] #  "diagonal"
    # Different initial volumes
    for waves in [3] # 1 2 
        speeds = zeros(Int, 4)
        if waves == 1
            speeds .= v_lam1_dia
        elseif waves == 2
            speeds .= v_lam2_dia
        elseif waves == 3
            speeds .= v_lam3_dia
        end
        for speed in [147] # speeds # 1 2 3
            pattern = "sine"
            println("Simulating moving substrate wettability with pattern $(pattern) and moving direction $(direction) and speed $(speed)")
            sys = Swalbe.SysConst(Lx=512, Ly=512, γ=0.01, δ=1.0, n=3, m=2, hmin=0.07, Tmax=5000000, tdump=5000)
            df_fluid = Dict()
            df_sub = Dict()
            θₚ = ones(sys.Lx,sys.Ly)
            # Substrate patterning
            if pattern == "sine" 
                for i in 1:sys.Lx, j in 1:sys.Ly
                    θₚ[i,j] = 1/9 + 1/18 * sin(2π*waves*(i-1)/sys.Lx) * sin(2π*waves*(j-1)/sys.Ly)
                end
            elseif pattern == "linear"
                θₚ = pyramidpattern(sys.Lx, sys.Ly, waves=waves)
            end
            # Make a cuarray with the substrate pattern
            θ_in = CUDA.adapt(CuArray, θₚ)
            if speed == 0
                # Actual simulation
                fluid, substrate = measure_substratewave(sys, "GPU", "blub", sub_speed=speed, θₛ=θ_in, dire=direction, dump=sys.tdump)
            else
                fluid, substrate = measure_substratewave(sys, "GPU", sub_speed=speed, θₛ=θ_in, dire=direction, dump=sys.tdump)
            end
            println("Writing measurements to Dict")
            # Filling the dataframes
            for t in 1:sys.Tmax÷sys.tdump
                df_fluid["h_$(t*sys.tdump)"] = fluid[t,:]
                df_sub["theta_$(t*sys.tdump)"] = substrate[t,:]
            end
            println("Saving Dict subdirection $direction subvel $speed and $(pattern) $waves to disk")
            bson("data/Moving_wettability/height_direc_$(direction)_sp_$(speed)_$(pattern)_$(waves)_tmax_$(sys.Tmax)_v2.bson", df_fluid)
            bson("data/Moving_wettability/theta_direc_$(direction)_sp_$(speed)_$(pattern)_$(waves)_tmax_$(sys.Tmax)_v2.bson", df_sub)
        
            CUDA.reclaim()
        end
    end
end

println("Script done, let's have a look at the data :)")
