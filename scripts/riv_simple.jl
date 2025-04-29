using DrWatson
@quickactivate :Swalbe
DataFrames, FileIO, Dates


# Fluid dynamics we need for the experiment
"""
    rivulet_run_sim

Simulation of a thin liquid rivulet with little noise added to the interface.

Checking the stability of a liquid rivulet under the influence of different forces.
What we consider here is the athermal or deterministic case as well as the fluctuating case.
The thermal fluctuations are considered using the fluctuating thin film equation. 

# Arguments
- `sys::Swalbe.SysConst`: Simulation relevant parameter, e.g. viscosity and slip
- `device::String`: Either "CPU" or "GPU"
- `R::AbstractFloat`: Radius in (x,y)-plane of the ring structure
- `rr::AbstractFloat`: Radius in (x,z)-plane
- `ϵ::AbstractFloat`: Amplititude of the initial interface undulations
- `dump::Int`: Dumping frequency for output creation  
- `fluid::AbstractArray`: Allocation of the output array (contains the results)
- `verbos::bool`: Switch to make the simulation write to console while running
- `ref2::bool`: Switch to use a different slip model

"""
function rivulet_run_sim(
    sys::Swalbe.SysConst;
    R = 150,
    rr = 100,
    ϵ = 0.01,
    dump = 1000,  
    fluid=zeros(sys.param.Tmax÷dump, sys.Lx*sys.Ly),
    verbos=true
)
    println("Running a simulation on rivulet stability\nThe rivulet is curved and resembles a torus\nThe substrate is uniform with theta = $(sys.param.θ)")
    state = Swalbe.Sys(sys, "CPU")
    # Set up initial condition
    
    h = Swalbe.torus(sys.Lx, sys.Ly, rr, R, sys.param.θ, (sys.Lx÷2, sys.Ly÷2), 0.12, noise=ϵ)
    state.height .= h
    
    println("Initial condition has been computed on the $(device)")
    Swalbe.equilibrium!(state, sys)
    state.ftemp .= state.fout
    println("Entering time loop")
    for t in 1:sys.param.Tmax
        if (t % sys.param.tdump == 0) || (t == 1)
            mass = 0.0
            mass = sum(state.height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
                println("Pressure $(sum(state.h∇px))\nSlip $(sum(state.slipx))")
            end
        end
        
        Swalbe.filmpressure!(state, sys)
        
        Swalbe.h∇p!(state)
        
        Swalbe.slippage_ring_riv!(state, sys)
        
        
        # Forces are 
        #   - pressure gradient 
        #   - substrate friction and slippage
        #   - thermal fluctuation (if sys.param.kbt > 0)
        state.Fx .= -state.h∇px .- state.slipx # .- state.kbtx
        state.Fy .= -state.h∇py .- state.slipy # .- state.kbty
        # println("Forces $(sum(state.Fx))\nP $(sum(state.pressure))\nh∇px $(sum(state.h∇px))\nSlip $(sum(state.slipx))")

        # New equilibrium
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        # New moments
        Swalbe.moments!(state)
        # Measurements, in this case only snapshots of simulational arrays
        Swalbe.snapshot!(fluid, state.height, t, dumping = dump)
    end
    return fluid
end

# Set up the simulation 
timeInterval = 1000

# Make a parameter sweep
for inC in [(80, 40, 1/9)] # 1/9, 1/6,  
    ang = inC[3]
    sys = Swalbe.SysConst(512, 512, Swalbe.Taumucs(Tmax=10000, δ=0.5, n=3, m=2, θ=ang, hmin=0.12)) # , hmin=0.3, hcrit=0.2
    # Run the simulation
    # Translate the initial conditions in seperate values
    outerRad = inC[1]
    innerRad = inC[2]
    # Run the simulation
    fluid = rivulet_run_sim(sys, R=outerRad, rr=innerRad, dump=timeInterval) #
    df_fluid = Dict()
    nSnapshots = sys.param.Tmax ÷ timeInterval
    for t in 1:nSnapshots
        # println("In saving loop at $(t) with $(size(fluid[t,:]))")
        df_fluid["h_$(t * timeInterval)"] = fluid[t,:]
    end
    println("Saving rivulet snapshots for R=$(outerRad) and r=$(innerRad) to disk")
    save_ang = round(Int,rad2deg(π*sys.param.θ))
    
    file_name = "data/Rivulets/test_riv.jld2"
   
    save(file_name, df_fluid)

    fluid .= 0.0
    df_fluid = Dict()
    println("Done with $(ang) $(outerRad) $(innerRad)")
end
