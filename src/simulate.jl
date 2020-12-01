# mutable struct dists{T} 
#     feq::Array{T, 3}
#     ftemp::Array{T, 3}
#     fout::Array{T, 3}
# end
"""
    run_flat(Sys::SysConst, device::String)

Performs a simulation of an flat interface without forces
"""
function run_flat(sys::SysConst, device::String)
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
        Swalbe.momen
end