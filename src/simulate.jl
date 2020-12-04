"""
    run_flat(Sys::SysConst, device::String)

Performs a simulation of an flat interface without forces
"""
function run_flat(sys::SysConst, device::String; verbos=true)
    println("Simulating a flat interface without driving forces (nothing should happen)")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    height .= 1.0
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % 100 == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(pressure, height, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, Fx, Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height
end

"""
    run_random(sys::SysConst, device::String)

Simulation of an random undulated interface
"""
function run_random(sys::SysConst, device::String; h₀=1.0, ϵ=0.01, verbos=true)
    println("Simulating a random undulated interface")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    Swalbe.randinterface!(height, h₀, ϵ)
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % 500 == 0
            mass = 0.0
            mass = sum(height)
            difference = maximum(height) - minimum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))\nAbsolute difference is $difference")
            end
        end
        Swalbe.filmpressure!(pressure, height, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height
end

# """
#     run_rayleightaylor(sys::SysConst, device::String)

# Simulation of an random undulated interface
# """
# function run_rayleightaylor(sys::SysConst, device::String; hinitial=ones(sys.Lx, sys.Ly), verbos=true)
#     println("Simulating the Rayleigh Taylor instability")
#     fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
#     # Swalbe.randinterface!(height, h₀, ϵ)
#     height .= hinitial
#     Swalbe.equilibrium!(feq, height, velx, vely, vsq)
#     ftemp .= feq
#     for t in 1:sys.Tmax
#         if t % 500 == 0
#             mass = 0.0
#             mass = sum(height)
#             difference = maximum(height) - minimum(height)
#             if verbos
#                 println("Time step $t mass is $(round(mass, digits=3))\nAbsolute difference is $difference")
#             end
#         end
#         Swalbe.filmpressure!(pressure, height, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
#         Swalbe.∇f!(h∇px, h∇py, pressure, height)
#         Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
#         Fx .= h∇px .+ slipx
#         Fy .= h∇py .+ slipy
#         Swalbe.equilibrium!(feq, height, velx, vely, vsq, sys.g)
#         Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
#         Swalbe.moments!(height, velx, vely, fout)
#     end
#     return height
# end

"""
    run_dropletrelax()

Simulates an out of equilibrium droplet
"""
function run_dropletrelax(sys::SysConst, device::String; radius=20, θ₀=1/6, center=(sys.Lx÷2, sys.Ly÷2), verbos=true)
    println("Simulating an out of equilibrium droplet")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    Swalbe.singledroplet(height, radius, θ₀, center)
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % 500 == 0
            mass = 0.0
            mass = sum(height)
            difference = maximum(height) - minimum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(pressure, height, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height
end

"""
    run_dropletpatterned()

Simulates an droplet on a patterned substrate
"""
function run_dropletpatterned(sys::SysConst, device::String; radius=20, θ₀=1/6, center=(sys.Lx÷2, sys.Ly÷2), θₛ=ones(sys.Lx, sys.Ly), verbos=true)
    println("Simulating a droplet on a patterned substrate")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    Swalbe.singledroplet(height, radius, θ₀, center)
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % 5000 == 0
            mass = 0.0
            mass = sum(height)
            maxU = maximum(abs.(velx))
            maxV = maximum(abs.(vely))
            if verbos
                println("Time step $t mass is $(round(mass, digits=3)) and max vel ($maxU $maxV)")
            end
        end
        Swalbe.filmpressure!(pressure, height, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height

end

"""
    run_dropletforced()

Simulates an droplet on a patterned substrate
"""
function run_dropletforced(sys::SysConst, device::String; radius=20, θ₀=1/6, center=(sys.Lx÷2, sys.Ly÷2), θₛ=fill(1/9, sys.Lx, sys.Ly), fx=0.0, fy=0.0, verbos=true)
    println("Simulating a sliding droplet")
    fout, ftemp, feq, height, velx, vely, vsq, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    Swalbe.singledroplet(height, radius, θ₀, center)
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % 1000 == 0
            mass = 0.0
            mass = sum(height)
            maxU = maximum(abs.(velx))
            maxV = maximum(abs.(vely))
            if verbos
                println("Time step $t mass is $(round(mass, digits=3)) and max vel ($maxU $maxV)")
            end
        end
        Swalbe.filmpressure!(pressure, height, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        # Here we a force that is like pull of an inclined plane
        Fx .= h∇px .+ slipx .+ fx .* height 
        Fy .= h∇py .+ slipy .+ fy .* height
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
    end
    return height, velx, velx, vely
    #= Works on the GPU =)
       sys = Swalbe.SysConst(Lx=512, Ly=512, Tmax=50000, δ=2.0)
       h = Swalbe.run_dropletforced(sys, "GPU", radius=50, θₛ=CUDA.fill(1/9, sys.Lx, sys.Ly), fx=-0.0001f0)
    =#
end
