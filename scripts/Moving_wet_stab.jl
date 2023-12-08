using DrWatson
@quickactivate :Swalbe
using CUDA, DataFrames, FileIO
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
    fluid=zeros(sys.param.Tmax÷dump, sys.Lx*sys.Ly),
    theta=zeros(sys.param.Tmax÷dump, sys.Lx*sys.Ly),
    dire = "x",
    verbos=true, 
    T=Float64
)
    state = Swalbe.Sys(sys, device)
    if device == "CPU"
        for i in 1:sys.Lx, j in 1:sys.Ly
            state.height[i,j] = h₀ + ϵ * sin(2π*wave_x*(i-1)/sys.Lx) * sin(2π*wave_y*(j-1)/sys.Ly)
        end
    elseif device == "GPU"
        h = zeros(size(state.height))
        for i in 1:sys.Lx, j in 1:sys.Ly
            h[i,j] = h₀ + ϵ * sin(2π*wave_x*(i-1)/sys.Lx) * sin(2π*wave_y*(j-1)/sys.Ly)
        end
        # theta = CUDA.zeros(Float64, sys.Lx, sys.Ly)
        state.height .= CUDA.adapt(CuArray, h)
    end
    Swalbe.equilibrium!(state, sys)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            minh = minimum(state.height)
            if minh < 0.05
                println("Rupture reached at $t")
                break
            end
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        
        Swalbe.filmpressure!(state, sys, θ=θₛ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        # Forces are the pressure gradient and the slippage due to substrate liquid boundary conditions
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        # New equilibrium
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        # New moments
        Swalbe.moments!(state)
        # Measurements, in this case only snapshots of simulational arrays
        Swalbe.snapshot!(fluid, state.height, t, dumping = dump)
        Swalbe.snapshot!(theta, θₛ, t, dumping = dump)
        move_substrate!(state.slipx, θₛ, t, sub_speed, direction=dire)
    end
    return fluid, theta
    if device == "GPU"
        CUDA.reclaim()
    end
end

function move_substrate!(θ, input, t, tmove; direction="diagonal")
    # No substrate velocity no moving
    if tmove == 0
        return
    # Finite substrate velocity
    else
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
end

direction = "diagonal"
# Different initial volumes
waves_nums = 20
theta_var = 1/18
run_on = "GPU"
speed = 0# 980
for waves_num in [20] # 4 5 6 7 8 9 10
    pattern = "sine"
    ang = 1/9
    TM = 600000 # 5000000
    println("Simulating moving substrate wettability with pattern $(pattern) and moving direction $(direction) and speed $(speed)")
    sys = Swalbe.SysConst(Lx=512, Ly=512, param=Swalbe.Taumucs(γ=0.01, δ=1.0, μ=1/12, n=3, m=2, hmin=0.07, Tmax=TM, tdump=1000))
    #sys = Swalbe.SysConst(Lx=512, Ly=512, γ=0.01, δ=1.0, n=3, m=2, hmin=0.07, Tmax=75000, tdump=500)
    df_fluid = Dict()
    df_sub = Dict()
    θₚ = ones(sys.Lx,sys.Ly)
    # Substrate patterning
    if pattern == "sine" 
        for i in 1:sys.Lx, j in 1:sys.Ly
            θₚ[i,j] = ang + theta_var * sin(2π*waves_num*(i-1)/sys.Lx) * sin(2π*waves_num*(j-1)/sys.Ly)
        end
    end
    # Make a cuarray with the substrate pattern
    θ_in = CUDA.adapt(CuArray, θₚ)
    fluid, substrate = measure_substratewave(sys, run_on, sub_speed=speed, θₛ=θ_in, dire=direction, dump=sys.param.tdump)
        
    println("Writing measurements to Dict")
    # Filling the dataframes
    for t in 1:sys.param.Tmax÷sys.param.tdump
        df_fluid["h_$(t*sys.param.tdump)"] = fluid[t,:]
        df_sub["theta_$(t*sys.param.tdump)"] = substrate[t,:]
    end
    println("Saving Dict subdirection $direction subvel $speed and $(pattern) $waves_num to disk")
    save_ang = Int(round(rad2deg(π*ang)))
    save_ang_del = Int(round(rad2deg(π*theta_var)))
    # save("data/Moving_wettability/height_direc_$(direction)_sp_$(speed)_$(pattern)_$(waves_num)_$(save_ang)_del_$(save_ang_del)_tmax_$(sys.param.Tmax)_v3.jld2", df_fluid)

    CUDA.reclaim()
end

println("Script done, let's have a look at the data :)")