"""
    time_loop(sys, state)

Time stepping procedure for the lattice Boltzmann state `state` given parameters `sys`
"""
function time_loop(sys::SysConst, state::LBM_state_2D; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    return state
end
function time_loop(sys::SysConst, state::LBM_state_2D, θ; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, θ=θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    return state
end

function time_loop(sys::SysConst, state::LBM_state_2D, Δh::Vector; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        push!(Δh, maximum(state.height) - minimum(state.height))
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    return state
end

function time_loop(sys::SysConst, state::LBM_state_2D, f::Function, measure::Vector; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        f(measure, state)
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    
    end
    return state, measure
end
# Time loop with snapshots saved
function time_loop(sys::SysConst, state::State, θ, data::Matrix; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, θ=θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= state.h∇py .+ state.slipy
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(data, state.height, t, dumping = sys.param.tdump)
        
    end
    return state
end

function time_loop(sys::SysConst_1D, state::State_1D; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    return state
end

function time_loop(sys::SysConst_1D, state::State_1D, θ; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, θ=θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    return state
end

function time_loop(sys::SysConst_1D, state::State_1D, Δh::Vector; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        push!(Δh, maximum(state.height) - minimum(state.height))
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    return state
end

function time_loop(sys::SysConst_1D, state::State_1D, f::Function, measure; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        f(measure, state)
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    
    end
    return state, measure
end

function time_loop(sys::SysConst_1D, state::State_1D, data::Matrix; verbose=false)
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(data, state.height, t, dumping = sys.param.tdump)
        
    end
    return data
end
"""
    run_flat(Sys::SysConst, device::String)

Performs a simulation of an flat interface without forces.

# Arguments

- `verbos :: Bool`: Enables consol output
- `T :: AbstractFloat`: Precision of output, default `Float64`

# Theory

Nothing at all should happen.
As the initial state is falt and no force is applied the fluid has no way to flow.
The equality
`` h(\\mathbf{x},0) = h(\\mathbf{x},\\infty), ``
should be satisfied for arbitrary many time steps.

# Examples
```julia
julia> using Swalbe, Test

julia> sys = Swalbe.SysConst(Lx=100, Ly=100, Tmax=5000);

julia> h = Swalbe.run_flat(sys, "CPU", verbos=false);

julia> @test all(h.height .== 1.0) # Check if all height values are identical to 1.0 (initial condition)
Test passed
``` 
"""
function run_flat(sys::SysConst, device::String; verbos=true)
    println("Simulating a flat interface without driving forces (nothing should happen) in two dimensions")
    state = Swalbe.Sys(sys, device)
    state.height .= 1.0
    time_loop(sys, state, verbose=verbos)

    return state.height
end

function run_flat(sys::SysConst_1D; verbos=true, T=Float64)
    println("Simulating a flat interface without driving forces (nothing should happen) in one dimension")
    state = Swalbe.Sys(sys)
    state.height .= 1.0
    time_loop(sys, state, verbose=verbos)

    return state.height
end

"""
    run_random(sys::SysConst, device::String)

Simulation of an random undulated interface

# Arguments
- `h₀ :: Float` : Average initial height
- `ϵ :: Float` : Amplitude of the flucutation
- `verbos :: Bool`: Enables consol output
- `T :: AbstractFloat`: Precision of output, default `Float64`

# Theory

Initial randomly perturbed fluid surface.
Unstable wavemodes should grow while wavemodes larger q₀ should be damped out.
Measuring this is on the TODO list.

# Examples
```julia
julia> using Swalbe, Test

julia> sys = Swalbe.SysConst(Lx=100, Ly=100, Tmax=5000);

julia> Swalbe.randinterface!(height, h₀, ϵ)

julia> h = Swalbe.run_random(sys, "CPU", h₀=10, ϵ=0.1, verbos=false);
```
"""
function run_random(sys::SysConst, device::String; h₀=1.0, ϵ=0.01, verbos=true)
    println("Simulating a random undulated interface in two dimensions")
    state = Swalbe.Sys(sys, device)
    Swalbe.randinterface!(state.height, h₀, ϵ)
    Swalbe.equilibrium!(state, sys)
    time_loop(sys, state, 1/9, verbose=verbos)
    
    return state.height
end
# 1D case
function run_random(sys::SysConst_1D; h₀=1.0, ϵ=0.01, verbos=true)
    println("Simulating a random undulated interface in one dimension")
    state = Swalbe.Sys(sys)
    Swalbe.randinterface!(state.height, h₀, ϵ)
    Swalbe.equilibrium!(state, sys)
    time_loop(sys, state, 1/9, verbose=verbos)
    return state.height
end


"""
    run_rayleightaylor(sys::SysConst, device::String)

Simulation of an random undulated interface and a gravitanional pull.

# Arguments
- `kx :: Int` : wavemode in x-direction, kx=18 -> 18 sine waves fitting into the domain
- `ky :: Int` : wavemode in y-direction
- `h₀ :: Float` : Average initial height
- `ϵ :: Float` : Amplitude of the flucutation
- `verbos :: Bool`: Enables consol output
- `T :: AbstractFloat`: Precision of output, default `Float64`

# Theory

Initial randomly perturbed fluid surface hanging from a substrate.
Here we have an interplay between gravity and surface tension.
The critical wavemode can be computed according to 
`` q_0 =  ``
Measuring this is on the TODO list.

# Examples
```julia
julia> using Swalbe, Test

julia> sys = Swalbe.SysConst(Lx=100, Ly=100, Tmax=5000);

julia> Swalbe.randinterface!(height, h₀, ϵ)

julia> h = Swalbe.run_random(sys, "CPU", h₀=10, ϵ=0.1, verbos=false);
```
"""
function run_rayleightaylor(sys::SysConst, device::String; kx=15, ky=18, h₀=1.0, ϵ=0.001, verbos=true, T=Float64)
    println("Simulating the Rayleigh Taylor instability in two dimensions")
    state = Swalbe.Sys(sys, device)
    for i in 1:sys.Lx, j in 1:sys.Ly
        state.height[i,j] = h₀ * (1 + ϵ * sin(2π*kx*i/(sys.Lx-1)) * sin(2π*ky*j/(sys.Ly-1)))
    end
    diff = []
    Swalbe.equilibrium!(state, sys)
    Swalbe.time_loop(sys, state, diff, verbose=verbos)
    return state.height, diff
end
# 1D case
function run_rayleightaylor(sys::SysConst_1D; k=15, h₀=1.0, ϵ=0.001, verbos=true, T=Float64)
    println("Simulating the Rayleigh Taylor instability in one dimension")
    state = Swalbe.Sys(sys)
    for i in 1:sys.L
        state.height[i] = h₀ * (1 + ϵ * sin(2π*k*i/(sys.L-1)))
    end
    diff = []
    Swalbe.equilibrium!(state, sys)
    Swalbe.time_loop(sys, state, diff, verbose=verbos)
    return state.height, diff
end

"""
    run_dropletrelax()

Simulates an out of equilibrium droplet
"""
function run_dropletrelax(sys::SysConst, device::String; radius=20, θ₀=1/6, center=(sys.Lx÷2, sys.Ly÷2), verbos=true)
    println("Simulating an out of equilibrium droplet in two dimensions")
    state = Swalbe.Sys(sys, device)
    Swalbe.singledroplet(state.height, radius, θ₀, center)
    Swalbe.equilibrium!(state, sys)
    area = []
    time_loop(sys, state, Swalbe.wetted!, area, verbose=verbos)
    return state.height, area
end
# 1D case
function run_dropletrelax(sys::SysConst_1D; radius=20, θ₀=1/6, center=(sys.L÷2), verbos=true )
    println("Simulating an out of equilibrium droplet in one dimensions")
    area = []
    state = Swalbe.Sys(sys)
    Swalbe.singledroplet(state.height, radius, θ₀, center) 
    Swalbe.equilibrium!(state, sys)
    time_loop(sys, state, Swalbe.wetted!, area, verbose=verbos)
    return state.height, area
end

"""
    run_dropletpatterned()

Simulates an droplet on a patterned substrate
"""
function run_dropletpatterned(sys::SysConst, device::String; radius=20, θ₀=1/6, center=(sys.Lx÷2, sys.Ly÷2), θₛ=fill(1/9, sys.Lx, sys.Ly), verbos=true)
    println("Simulating a droplet on a patterned substrate in two dimensions")
    state = Swalbe.Sys(sys, device)
    Swalbe.singledroplet(state.height, radius, θ₀, center)
    Swalbe.equilibrium!(state, sys)
    time_loop(sys, state, θₛ, verbose=verbos)
    
    return state.height

end
# 1D case
function run_dropletpatterned(sys::SysConst_1D; radius=20, θ₀=1/6, center=sys.L÷2, θₛ=fill(1/2, sys.L), verbos=true)
    println("Simulating a droplet on a patterned substrate in one dimension")
    state = Swalbe.Sys(sys)
    Swalbe.singledroplet(state.height, radius, θ₀, center)
    Swalbe.equilibrium!(state, sys)
    time_loop(sys, state, θₛ, verbose=verbos)
    return state.height
end

"""
    run_dropletforced()

Simulates an droplet on a patterned substrate
"""
function run_dropletforced(sys::SysConst, device::String; radius=20, θ₀=1/6, center=(sys.Lx÷2, sys.Ly÷2), fx=0.0, fy=0.0, verbos=true)
    bodyforce = zeros(2)
    bodyforce[1] = fx
    bodyforce[2] = fy
    println("Simulating a sliding droplet in two dimensions")
    state = Swalbe.Sys(sys, device)
    Swalbe.singledroplet(state.height, radius, θ₀, center)
    Swalbe.equilibrium!(state, sys)
    println("Starting the lattice Boltzmann time loop")
    time_loop(sys, state, Swalbe.inclination!, bodyforce, verbose=verbos)

    return state.height, state.velx, state.vely
    
end
# 1D case
function run_dropletforced(sys::SysConst_1D; radius=20, θ₀=1/6, center=(sys.L÷2), θₛ=fill(1/9, sys.L), f=0.0, verbos=true)
    println("Simulating a sliding droplet in one dimension")
    state = Swalbe.Sys(sys)
    Swalbe.singledroplet(state.height, radius, θ₀, center)
    Swalbe.equilibrium!(state, sys)
    println("Starting the lattice Boltzmann time loop")
    time_loop(sys, state, Swalbe.inclination!, f, verbose=verbos)
    return state.height, state.vel

end