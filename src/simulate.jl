"""
    run_flat(Sys::SysConst, device::String)

Performs a simulation of an flat interface without forces
"""
function run_flat(sys::SysConst, device::String)
    println("Simulating a flat interface without driving forces (nothing should happen)")
    fout, ftemp, feq, height, velx, vely, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    height .= 1.0
    vsq = zeros(size(height))
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    ftemp .= feq
    for t in 1:sys.Tmax
        if t % 100 == 0
            mass = 0.0
            mass = sum(height)
            println("Time step $t mass is $(round(mass, digits=3))")
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
    println("Simulating a flat interface without driving forces (nothing should happen)")
    fout, ftemp, feq, height, velx, vely, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false)
    Swalbe.randinterface!(height, h₀, ϵ)
    vsq = zeros(size(height))
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
