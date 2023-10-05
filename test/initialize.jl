@testset "Allocations" begin
    T = Float64
    sys = Swalbe.SysConst{T}(Lx=25, Ly=26, param=Swalbe.Taumucs())
    @test sys.Lx == 25
    @test sys.Ly == 26
    @test sys.param.Tmax == 1000
    @test sys.param.tdump == 100

    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py, fthermalx, fthermaly = Swalbe.Sys(sys, "CPU", true, T)
    statedict = Dict("fo"=>fout, 
                     "ft"=>ftemp, 
                     "fe"=>feq, 
                     "h"=>height, 
                     "vx"=>velx, 
                     "vy"=>vely, 
                     "vv"=>vsq, 
                     "p"=>pressure, 
                     "grad"=>dgrad, 
                     "fx"=>Fx, 
                     "fy"=>Fy, 
                     "sx"=>slipx, 
                     "sy"=>slipy, 
                     "px"=>h∇px, 
                     "py"=>h∇py, 
                     "tx"=>fthermalx, 
                     "ty"=>fthermaly)

    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, "CPU", false, T)
    statedict2 = Dict("fo"=>fout, 
                      "ft"=>ftemp, 
                      "fe"=>feq, 
                      "h"=>height, 
                      "vx"=>velx, 
                      "vy"=>vely, 
                      "vv"=>vsq, 
                      "p"=>pressure, 
                      "grad"=>dgrad, 
                      "fx"=>Fx, 
                      "fy"=>Fy, 
                      "sx"=>slipx, 
                      "sy"=>slipy, 
                      "px"=>h∇px, 
                      "py"=>h∇py)

    dyn = Swalbe.Sys(sys, "CPU")
    @test isa(dyn, Swalbe.State)
    dyndict = Dict("fo"=>dyn.fout, 
                   "ft"=>dyn.ftemp, 
                   "fe"=>dyn.feq, 
                   "h"=>dyn.height, 
                   "vx"=>dyn.velx, 
                   "vy"=>dyn.vely, 
                   "vv"=>dyn.vsq, 
                   "p"=>dyn.pressure, 
                   "grad"=>dyn.dgrad, 
                   "fx"=>dyn.Fx, 
                   "fy"=>dyn.Fy, 
                   "sx"=>dyn.slipx, 
                   "sy"=>dyn.slipy, 
                   "px"=>dyn.h∇px, 
                   "py"=>dyn.h∇py)
    # Struct should do what I expect, hopefully
    dynth = Swalbe.Sys(sys, "CPU", kind="thermal")
    @test isa(dynth, Swalbe.State_thermal)
    dynthdict = Dict("fo"=>dynth.basestate.fout, 
                   "ft"=>dynth.basestate.ftemp, 
                   "fe"=>dynth.basestate.feq, 
                   "h"=>dynth.basestate.height, 
                   "vx"=>dynth.basestate.velx, 
                   "vy"=>dynth.basestate.vely, 
                   "vv"=>dynth.basestate.vsq, 
                   "p"=>dynth.basestate.pressure, 
                   "grad"=>dynth.basestate.dgrad, 
                   "fx"=>dynth.basestate.Fx, 
                   "fy"=>dynth.basestate.Fy, 
                   "sx"=>dynth.basestate.slipx, 
                   "sy"=>dynth.basestate.slipy, 
                   "px"=>dynth.basestate.h∇px, 
                   "py"=>dynth.basestate.h∇py,
                   "tx"=>dynth.kbtx,
                   "ty"=>dynth.kbty)

    for i in [statedict, statedict2, dyndict, dynthdict]
        @test isa(i["fo"], Array{Float64, 3})
        @test isa(i["ft"], Array{Float64, 3})
        @test isa(i["fe"], Array{Float64, 3})
        @test isa(i["h"], Array{Float64, 2})
        @test isa(i["vx"], Array{Float64, 2})
        @test isa(i["vy"], Array{Float64, 2})
        @test isa(i["vv"], Array{Float64, 2})
        @test isa(i["p"], Array{Float64, 2})
        @test isa(i["grad"], Array{Float64, 3})
        @test isa(i["fx"], Array{Float64, 2})
        @test isa(i["fy"], Array{Float64, 2})
        @test isa(i["sx"], Array{Float64, 2})
        @test isa(i["sy"], Array{Float64, 2})
        @test isa(i["px"], Array{Float64, 2})
        @test isa(i["py"], Array{Float64, 2})
        if (i == statedict) || (i == dynthdict)
            @test isa(i["tx"], Array{Float64, 2})
            @test isa(i["ty"], Array{Float64, 2})
        end
    end
    
    # 1D case
    sys1d = Swalbe.SysConst_1D{T}(L=25, param=Swalbe.Taumucs())
    @test sys1d.L == 25
    @test sys1d.param.Tmax == 1000
    @test sys1d.param.tdump == 100
    dy1 = Swalbe.Sys(sys1d)
    @test typeof(dy1) == Swalbe.State_1D{Float64}
    @test typeof(dy1) <: Swalbe.LBM_state_1D
    dydict = Dict("fo"=>dy1.fout, 
                  "ft"=>dy1.ftemp, 
                  "fe"=>dy1.feq, 
                  "h"=>dy1.height, 
                  "v"=>dy1.vel,
                  "p"=>dy1.pressure, 
                  "grad"=>dy1.dgrad, 
                  "f"=>dy1.F, 
                  "s"=>dy1.slip,
                  "pg"=>dy1.h∇p)
    dy2 = Swalbe.Sys(sys1d, kind="thermal")
    @test typeof(dy2) == Swalbe.State_thermal_1D{Float64}
    @test typeof(dy2) <: Swalbe.LBM_state_1D
    dydict2 = Dict("fo"=>dy2.basestate.fout, 
                  "ft"=>dy2.basestate.ftemp, 
                  "fe"=>dy2.basestate.feq, 
                  "h"=>dy2.basestate.height, 
                  "v"=>dy2.basestate.vel,
                  "p"=>dy2.basestate.pressure, 
                  "grad"=>dy2.basestate.dgrad, 
                  "f"=>dy2.basestate.F, 
                  "s"=>dy2.basestate.slip,
                  "pg"=>dy2.basestate.h∇p,
                  "t"=>dy2.kbt)
    dy3 = Swalbe.Sys(sys1d, kind="gamma")
    @test typeof(dy3) == Swalbe.State_gamma_1D{Float64}
    @test typeof(dy3) <: Swalbe.LBM_state_1D
    dydict3 = Dict("fo"=>dy3.basestate.fout, 
                  "ft"=>dy3.basestate.ftemp, 
                  "fe"=>dy3.basestate.feq, 
                  "h"=>dy3.basestate.height, 
                  "v"=>dy3.basestate.vel,
                  "p"=>dy3.basestate.pressure, 
                  "grad"=>dy3.basestate.dgrad, 
                  "f"=>dy3.basestate.F, 
                  "s"=>dy3.basestate.slip,
                  "pg"=>dy3.basestate.h∇p,
                  "ga"=>dy3.γ,
                  "dga"=>dy3.∇γ)
    # Struct should do what I expect, hopefully
    for i in [dydict, dydict2, dydict3]
        @test isa(i["fo"], Matrix{Float64})
        @test isa(i["ft"], Matrix{Float64})
        @test isa(i["fe"], Matrix{Float64})
        @test isa(i["h"], Vector{Float64})
        @test isa(i["v"], Vector{Float64})
        @test isa(i["p"], Vector{Float64})
        @test isa(i["grad"], Matrix{Float64})
        @test isa(i["f"], Vector{Float64})
        @test isa(i["s"], Vector{Float64})
        @test isa(i["pg"], Vector{Float64})
        if i == dydict2
            @test isa(i["t"], Vector{Float64})
        end
        if i == dydict3
            @test isa(i["ga"], Vector{Float64})
            @test isa(i["dga"], Vector{Float64})
        end
    end
end
