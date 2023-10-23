using Swalbe
using Plots, CUDA, DataFrames, FileIO
# CUDA.device!(1)

# Fluid dynamics we need for the experiment
"""
    torus_stability

Simulates a torus like fluid structure.
"""
function torus_stability(
    sys::Swalbe.SysConst, 
    device::String; 
    R = 150,
    rr = 100,
    h₀ = 2.0, 
    ϵ = 0.01,
    dump = 1000,  
    fluid=zeros(sys.param.Tmax÷dump, sys.Lx*sys.Ly),
    verbos=true, 
    T=Float64
)
    println("Running a simulation on rivulet stability\nThe rivulet is curved and resembles a torus")
    state = Swalbe.Sys(sys, device)
    # Set up initial condition
    h = Swalbe.torus(sys.Lx, sys.Ly, rr, R, sys.param.θ, (sys.Lx÷2, sys.Ly÷2))
    # Push it to the desired device
    if device == "CPU"
        state.height .= h
    elseif device == "GPU"
        CUDA.copyto!(state.height, h)
    end
    println("Initial condition has been computed on the $(device)")
    Swalbe.equilibrium!(state, sys)
    state.ftemp .= state.fout
    println("Entering time loop")
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        
        Swalbe.filmpressure!(state, sys)
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
    end
    return fluid
    CUDA.reclaim()
end
