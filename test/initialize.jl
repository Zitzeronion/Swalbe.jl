@testset "Allocations" begin
    T = Float64
    sys = Swalbe.SysConst{T}(Lx=25, Ly=26)
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py, fthermalx, fthermaly = Swalbe.Sys(sys, "CPU", true, T)
    # Struct should do what I expect, hopefully
    @test sys.Lx == 25
    @test sys.Ly == 26
    @test sys.Tmax == 1000
    @test sys.tdump == 100
    @test isa(fout, Array{Float64, 3})
    @test isa(ftemp, Array{Float64, 3})
    @test isa(feq, Array{Float64, 3})
    @test isa(height, Array{Float64, 2})
    @test isa(velx, Array{Float64, 2})
    @test isa(vely, Array{Float64, 2})
    @test isa(vsq, Array{Float64, 2})
    @test isa(pressure, Array{Float64, 2})
    @test isa(dgrad, Array{Float64, 3})
    @test isa(Fx, Array{Float64, 2})
    @test isa(Fy, Array{Float64, 2})
    @test isa(slipx, Array{Float64, 2})
    @test isa(slipy, Array{Float64, 2})
    @test isa(h∇px, Array{Float64, 2})
    @test isa(h∇py, Array{Float64, 2})
    @test isa(fthermalx, Array{Float64, 2})
    @test isa(fthermaly, Array{Float64, 2})
    
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py, = Swalbe.Sys(sys, "CPU", false, T)
    # Struct should do what I expect, hopefully
    @test sys.Lx == 25
    @test sys.Ly == 26
    @test isa(fout, Array{Float64, 3})
    @test isa(ftemp, Array{Float64, 3})
    @test isa(feq, Array{Float64, 3})
    @test isa(height, Array{Float64, 2})
    @test isa(velx, Array{Float64, 2})
    @test isa(vely, Array{Float64, 2})
    @test isa(vsq, Array{Float64, 2})
    @test isa(pressure, Array{Float64, 2})
    @test isa(dgrad, Array{Float64, 3})
    @test isa(Fx, Array{Float64, 2})
    @test isa(Fy, Array{Float64, 2})
    @test isa(slipx, Array{Float64, 2})
    @test isa(slipy, Array{Float64, 2})
    @test isa(h∇px, Array{Float64, 2})
    @test isa(h∇py, Array{Float64, 2})
    
    dyn = Swalbe.Sys(sys, "CPU")
    # Struct should do what I expect, hopefully
    @test sys.Lx == 25
    @test sys.Ly == 26
    @test isa(dyn.fout, Array{Float64, 3})
    @test isa(dyn.ftemp, Array{Float64, 3})
    @test isa(dyn.feq, Array{Float64, 3})
    @test isa(dyn.height, Array{Float64, 2})
    @test isa(dyn.velx, Array{Float64, 2})
    @test isa(dyn.vely, Array{Float64, 2})
    @test isa(dyn.vsq, Array{Float64, 2})
    @test isa(dyn.pressure, Array{Float64, 2})
    @test isa(dyn.dgrad, Array{Float64, 3})
    @test isa(dyn.Fx, Array{Float64, 2})
    @test isa(dyn.Fy, Array{Float64, 2})
    @test isa(dyn.slipx, Array{Float64, 2})
    @test isa(dyn.slipy, Array{Float64, 2})
    @test isa(dyn.h∇px, Array{Float64, 2})
    @test isa(dyn.h∇py, Array{Float64, 2})
    
    # 1D case
    sys1d = Swalbe.SysConst_1D{T}(L=25)
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p, fthermal = Swalbe.Sys(sys1d, true, T)
    # Struct should do what I expect, hopefully
    @test sys1d.L == 25
    @test sys1d.Tmax == 1000
    @test sys1d.tdump == 100
    @test isa(fout, Matrix{Float64})
    @test isa(ftemp, Matrix{Float64})
    @test isa(feq, Matrix{Float64})
    @test isa(height, Vector{Float64})
    @test isa(vel, Vector{Float64})
    @test isa(pressure, Vector{Float64})
    @test isa(dgrad, Matrix{Float64})
    @test isa(F, Vector{Float64})
    @test isa(slip, Vector{Float64})
    @test isa(h∇p, Vector{Float64})
    @test isa(fthermal, Vector{Float64})
    
    dyn = Swalbe.Sys(sys1d)
    # Struct should do what I expect, hopefully
    @test sys1d.L == 25
    @test sys1d.Tmax == 1000
    @test sys1d.tdump == 100
    @test isa(dyn.fout, Matrix{Float64})
    @test isa(dyn.ftemp, Matrix{Float64})
    @test isa(dyn.feq, Matrix{Float64})
    @test isa(dyn.height, Vector{Float64})
    @test isa(dyn.vel, Vector{Float64})
    @test isa(dyn.pressure, Vector{Float64})
    @test isa(dyn.dgrad, Matrix{Float64})
    @test isa(dyn.F, Vector{Float64})
    @test isa(dyn.slip, Vector{Float64})
    @test isa(dyn.h∇p, Vector{Float64})
    
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys1d, false, T)
    # Struct should do what I expect, hopefully
    @test sys1d.L == 25
    @test sys1d.Tmax == 1000
    @test sys1d.tdump == 100
    @test isa(fout, Matrix{Float64})
    @test isa(ftemp, Matrix{Float64})
    @test isa(feq, Matrix{Float64})
    @test isa(height, Vector{Float64})
    @test isa(vel, Vector{Float64})
    @test isa(pressure, Vector{Float64})
    @test isa(dgrad, Matrix{Float64})
    @test isa(F, Vector{Float64})
    @test isa(slip, Vector{Float64})
    @test isa(h∇p, Vector{Float64})
end
