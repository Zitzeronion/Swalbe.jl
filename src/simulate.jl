# """
#     run()

# Combine testcase into a single a run() function.
# WIP
# """
# function run(def::String, sys::SysConst, device::String; verbos=true)
#     dyn = Swalbe.Sys(sys, device) 
#     if def == "Rayleigh-Taylor"

#     if def == "Relax-droplet"

#     if def == "Moving-droplet"

#     end
# end

"""
    time_loop(sys, state, θ)

Time stepping procedure for the lattice Boltzmann state `state` given parameters `sys`
"""
function time_loop(sys::SysConst, state::State, θ; verbose=false)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
    end
    return state
end

function time_loop(sys::SysConst, state::State, θ, Δh::Vector; verbose=false)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Δh[t] = maximum(h) - minimum(h)
        Swalbe.filmpressure!(state, sys, θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= -state.h∇py .- state.slipy
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
    end
    return state
end
# Time loop with snapshots saved
function time_loop(sys::SysConst, state::State, θ, data::Matrix; verbose=false)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.Fx .= -state.h∇px .- state.slipx
        state.Fy .= state.h∇py .+ state.slipy
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(data, state.height, t, dumping = sys.tdump)
        
    end
    return state
end

function time_loop(sys::SysConst_1D, state::State_1D, θ; verbose=false)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, θ)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
    end
    return state
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
    time_loop(sys, state, 1/9, verbose=verbos)

    return state.height
end

function run_flat(sys::SysConst_1D; verbos=true, T=Float64)
    println("Simulating a flat interface without driving forces (nothing should happen) in one dimension")
    state = Swalbe.Sys(sys)
    state.height .= 1.0
    time_loop(sys, state, 1/9, verbose=verbos)

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
    Swalbe.equilibrium!(state)
    time_loop(sys, state, 1/9, verbose=verbos)
    
    return state.height
end
# 1D case
function run_random(sys::SysConst_1D; h₀=1.0, ϵ=0.01, verbos=true)
    println("Simulating a random undulated interface in one dimension")
    state = Swalbe.Sys(sys)
    Swalbe.randinterface!(state.height, h₀, ϵ)
    Swalbe.equilibrium!(state)
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
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    for i in 1:sys.Lx, j in 1:sys.Ly
        height[i,j] = h₀ * (1 + ϵ * sin(2π*kx*i/(sys.Lx-1)) * sin(2π*ky*j/(sys.Ly-1)))
    end
    diff = []
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax 
        push!(diff, maximum(height) - minimum(height))
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))\nAbsolute difference is $difference")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq, sys.g)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height, diff
end
# 1D case
function run_rayleightaylor(sys::SysConst_1D; k=15, h₀=1.0, ϵ=0.001, verbos=true, T=Float64)
    println("Simulating the Rayleigh Taylor instability in one dimension")
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys, false, T)
    for i in 1:sys.L
        height[i] = h₀ * (1 + ϵ * sin(2π*k*i/(sys.L-1)))
    end
    diff = []
    Swalbe.equilibrium!(feq, height, vel)
    ftemp .= feq
    for t in 1:sys.Tmax
        push!(diff, maximum(height) - minimum(height))
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))\nAbsolute difference is $difference")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        F .= h∇p .+ slip
        Swalbe.equilibrium!(feq, height, vel, sys.g)
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        Swalbe.moments!(height, vel, fout)
    end
    return height, diff
end

"""
    run_dropletrelax()

Simulates an out of equilibrium droplet
"""
function run_dropletrelax(
    sys::SysConst, 
    device::String; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.Lx÷2, sys.Ly÷2), 
    verbos=true, 
    T=Float64
)
    println("Simulating an out of equilibrium droplet in two dimensions")
    area = []
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    if device == "CPU"
        Swalbe.singledroplet(height, radius, θ₀, center)
    elseif device == "GPU"
        h = zeros(size(height))
        Swalbe.singledroplet(h, radius, θ₀, center)
        height = CUDA.adapt(CuArray, h)
    end
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        push!(area, length(findall(height .> 0.055)))
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height, area
end
# 1D case
function run_dropletrelax(
    sys::SysConst_1D;
    radius=20, 
    θ₀=1/6, 
    center=(sys.L÷2), 
    verbos=true, 
    T=Float64
)
    println("Simulating an out of equilibrium droplet in one dimensions")
    area = []
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys, false, T)
    
    Swalbe.singledroplet(height, radius, θ₀, center)
    
    Swalbe.equilibrium!(feq, height, vel)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        push!(area, length(findall(height .> 0.055)))
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        F .= h∇p .+ slip
        Swalbe.equilibrium!(feq, height, vel)
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        Swalbe.moments!(height, vel, fout)
    end
    return height, area
end

"""
    run_dropletpatterned()

Simulates an droplet on a patterned substrate
"""
function run_dropletpatterned(
    sys::SysConst, 
    device::String; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.Lx÷2, sys.Ly÷2), 
    θₛ=ones(sys.Lx, sys.Ly), 
    verbos=true, 
    T=Float64
)
    println("Simulating a droplet on a patterned substrate in two dimensions")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    if device == "CPU"
        Swalbe.singledroplet(height, radius, θ₀, center)
    elseif device == "GPU"
        h = zeros(size(height))
        Swalbe.singledroplet(h, radius, θ₀, center)
        height = CUDA.adapt(CuArray, h)
    end
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            maxU = maximum(abs.(velx))
            maxV = maximum(abs.(vely))
            if verbos
                println("Time step $t mass is $(round(mass, digits=3)) and max vel ($maxU $maxV)")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height

end
# 1D case
function run_dropletpatterned(
    sys::SysConst_1D; 
    radius=20, 
    θ₀=1/6, 
    center=sys.L÷2, 
    θₛ=ones(sys.L), 
    verbos=true, 
    T=Float64
)
    println("Simulating a droplet on a patterned substrate in one dimension")
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys, false, T)
    
    Swalbe.singledroplet(height, radius, θ₀, center)
        
    Swalbe.equilibrium!(feq, height, vel)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            maxU = maximum(abs.(vel))
            if verbos
                println("Time step $t mass is $(round(mass, digits=3)) and max vel $maxU")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        F .= h∇p .+ slip
        Swalbe.equilibrium!(feq, height, vel)
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        Swalbe.moments!(height, vel, fout)
    end
    return height

end

"""
    run_dropletforced()

Simulates an droplet on a patterned substrate
"""
function run_dropletforced(
    sys::SysConst, 
    device::String; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.Lx÷2, sys.Ly÷2), 
    θₛ=fill(1/9, sys.Lx, sys.Ly), 
    fx=0.0, 
    fy=0.0, 
    verbos=true, 
    T=Float64
)
    println("Simulating a sliding droplet in two dimensions")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    if device == "CPU"
        Swalbe.singledroplet(height, radius, θ₀, center)
    elseif device == "GPU"
        h = zeros(size(height))
        Swalbe.singledroplet(h, radius, θ₀, center)
        height = CUDA.adapt(CuArray, h)
    end
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            maxU = maximum(abs.(velx))
            maxV = maximum(abs.(vely))
            if verbos
                println("Time step $t mass is $(round(mass, digits=3)) and max vel ($maxU $maxV)")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        # Here we a force that is like pull of an inclined plane
        Fx .= h∇px .+ slipx .+ fx .* height 
        Fy .= h∇py .+ slipy .+ fy .* height
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height, velx, vely
    #= Works on the GPU =)
       sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
       h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
    =#
end
# 1D case
function run_dropletforced(
    sys::SysConst_1D; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.L÷2), 
    θₛ=fill(1/9, sys.L), 
    f=0.0, 
    verbos=true, 
    T=Float64
)
    println("Simulating a sliding droplet in one dimension")
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys, false, T)
    
    Swalbe.singledroplet(height, radius, θ₀, center)
    
    Swalbe.equilibrium!(feq, height, vel)
    ftemp .= feq
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            maxU = maximum(abs.(vel))
            if verbos
                println("Time step $t mass is $(round(mass, digits=3)) and max vel ($maxU)")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        # Here we a force that is like pull of an inclined plane
        F .= h∇p .+ slip .+ f .* height 
        Swalbe.equilibrium!(feq, height, vel)
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        Swalbe.moments!(height, vel, fout)
    end
    return height, vel
    #= Works on the GPU =)
       sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
       h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
    =#
end

# """
#     run_fluctuating_thin_film()


# """
# function run_fluctuating_thin_film(
#     sys::SysConst, 
#     device::String;
#     h₀=1,
#     ϵ=0.001, 
#     θ₀=1/9,  
#     verbos=true, 
#     T=Float64
# )
#     println("Simulating a thin film with thermal fluctuations")
#     fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py, tx, ty = Swalbe.Sys(sys, device, true, T)
#     if device == "CPU"
#         Swalbe.randinterface!(height, h₀, ϵ)
#     elseif device == "GPU"
#         h = zeros(size(height))
#         Swalbe.randinterface!(h, h₀, ϵ)
#         height = CUDA.adapt(CuArray, h)
#     end
#     Swalbe.equilibrium!(feq, height, velx, vely, vsq)
#     ftemp .= feq
#     println("Starting the lattice Boltzmann time loop")
#     for t in 1:sys.Tmax
#         if t % sys.tdump == 0
#             mass = 0.0
#             mass = sum(height)
#             if verbos
#                 println("Time step $t mass is $(round(mass, digits=3))")
#             end
#         end
#         Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit)
#         Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
#         Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
#         Swalbe.thermal!(tx, ty, height, sys.kbt, sys.μ, sys.δ)
#         # Here we a force that is like pull of an inclined plane
#         Fx .= h∇px .+ slipx .+ tx 
#         Fy .= h∇py .+ slipy .+ ty
#         Swalbe.equilibrium!(feq, height, velx, vely, vsq)
#         Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
#         Swalbe.moments!(height, velx, vely, fout)
#     end
#     return height, velx, vely
#     #= Works on the GPU =)
#        sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
#        h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
#     =#
# end

# function run_fluctuating_thin_film(
#     sys::SysConst_1D;
#     h₀=1,
#     ϵ=0.001, 
#     θ₀=1/9,  
#     verbos=true, 
#     T=Float64
# )
#     println("Simulating a thin film with thermal fluctuations")
#     fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p, fluc = Swalbe.Sys(sys, true, T)
#     Swalbe.randinterface!(height, h₀, ϵ)
    
#     Swalbe.equilibrium!(feq, height, vel)
#     ftemp .= feq
#     println("Starting the lattice Boltzmann time loop")
#     for t in 1:sys.Tmax
#         if t % sys.tdump == 0
#             mass = 0.0
#             mass = sum(height)
#             if verbos
#                 println("Time step $t mass is $(round(mass, digits=3))")
#             end
#         end
#         Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit)
#         Swalbe.∇f!(h∇p, pressure, dgrad, height)
#         Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
#         Swalbe.thermal!(fluc, height, sys.kbt, sys.μ, sys.δ)
#         # Here we a force that is like pull of an inclined plane
#         F .= h∇p .+ slip .+ fluc 
#         Swalbe.equilibrium!(feq, height, vel)
#         Swalbe.BGKandStream!(fout, feq, ftemp, -F)
#         Swalbe.moments!(height, vel, fout)
#     end
#     return height, vel
#     #= Works on the GPU =)
#        sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
#        h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
#     =#
# end

# function run_no_fluctuating_thin_film(
#     sys::SysConst, 
#     device::String;
#     h₀=1,
#     ϵ=0.001, 
#     θ₀=1/9,  
#     verbos=true, 
#     T=Float64
# )
#     println("Simulating a thin film with thermal fluctuations")
#     fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
#     if device == "CPU"
#         Swalbe.randinterface!(height, h₀, ϵ)
#     elseif device == "GPU"
#         h = zeros(size(height))
#         Swalbe.randinterface!(h, h₀, ϵ)
#         height = CUDA.adapt(CuArray, h)
#     end
#     Swalbe.equilibrium!(feq, height, velx, vely, vsq)
#     ftemp .= feq
#     println("Starting the lattice Boltzmann time loop")
#     for t in 1:sys.Tmax
#         if t % sys.tdump == 0
#             mass = 0.0
#             mass = sum(height)
#             if verbos
#                 println("Time step $t mass is $(round(mass, digits=3))")
#             end
#         end
#         Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit)
#         Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
#         Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
#         # Swalbe.thermal!(tx, ty, height, sys.kbt, sys.μ, sys.δ)
#         # Here we a force that is like pull of an inclined plane
#         Fx .= h∇px .+ slipx # .+ tx 
#         Fy .= h∇py .+ slipy # .+ ty
#         Swalbe.equilibrium!(feq, height, velx, vely, vsq)
#         Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
#         Swalbe.moments!(height, velx, vely, fout)
#     end
#     return height, velx, vely
#     #= Works on the GPU =)
#        sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
#        h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
#     =#
# end

# function run_no_fluctuating_thin_film(
#     sys::SysConst_1D;
#     h₀=1,
#     ϵ=0.001, 
#     θ₀=1/9,  
#     verbos=true, 
#     T=Float64
# )
#     println("Simulating a thin film with thermal fluctuations")
#     fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys, false, T)
    
#     Swalbe.randinterface!(height, h₀, ϵ)
    
#     Swalbe.equilibrium!(feq, height, vel)
#     ftemp .= feq
#     println("Starting the lattice Boltzmann time loop")
#     for t in 1:sys.Tmax
#         if t % sys.tdump == 0
#             mass = 0.0
#             mass = sum(height)
#             if verbos
#                 println("Time step $t mass is $(round(mass, digits=3))")
#             end
#         end
#         Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit)
#         Swalbe.∇f!(h∇p, pressure, dgrad, height)
#         Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        
#         # Here we a force that is like pull of an inclined plane
#         F .= h∇p .+ slip # .+ tx 
        
#         Swalbe.equilibrium!(feq, height, vel)
#         Swalbe.BGKandStream!(fout, feq, ftemp, -F)
#         Swalbe.moments!(height, vel, fout)
#     end
#     return height, vel
#     #= Works on the GPU =)
#        sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
#        h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
#     =#
# end

# function run_active_thin_film(
#     sys::SysConst_1D;
#     h₀=1,
#     ϵ=0.001, 
#     θ₀=1/9,  
#     verbos=true, 
#     rho = zeros(sys.L),
#     lap_rho = zeros(sys.L), 
#     grad_rho = zeros(sys.L), 
#     grad_h = zeros(sys.L), 
#     lap_h = zeros(sys.L),
#     Gam = 0.5,
#     mobl = 1.0, 
#     difu = 0.1,
#     T=Float64
# )
#     println("Simulating a active thin film")
#     fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p, fluc = Swalbe.Sys(sys, true, T)
#     Swalbe.randinterface!(height, h₀, ϵ)
#     # Swalbe.sinewave1D!(height, h₀, 1, ϵ, 1)
#     rho .= 0.1
#     Swalbe.equilibrium!(feq, height, vel)
#     ftemp .= feq
#     println("Starting the lattice Boltzmann time loop")
#     for t in 1:sys.Tmax
#         if t % sys.tdump == 0
#             mass = 0.0
#             mass = sum(height)
#             if verbos
#                 println("Time step $t mass is $(round(mass, digits=3))")
#             end
#         end
#         Swalbe.update_rho(fluc, rho, dgrad, height, lap_rho, grad_rho, grad_h, lap_h, M=mobl, D=difu)
#         # Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit)
#         Swalbe.filmpressure!(pressure, height, dgrad, rho, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit, Gamma=Gam)
#         Swalbe.∇f!(h∇p, pressure, dgrad, height)
#         Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
#         # HereSwalbe.thermal!(fluc, height, sys.kbt, sys.μ, sys.δ)
#         # Here we a force that is like pull of an inclined plane
#         F .= h∇p .+ slip # .+ fluc 
#         Swalbe.equilibrium!(feq, height, vel)
#         Swalbe.BGKandStream!(fout, feq, ftemp, -F)
#         Swalbe.moments!(height, vel, fout)
#     end
    

#     return height, vel, rho
#
# end