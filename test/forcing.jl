@testset "Forcings" begin
    fx = zeros(5,5)
    fy = zeros(5,5)
    f1 = zeros(30)
    # For structs
    sys = Swalbe.SysConst(Lx=5, Ly=5)
    sys_act = Swalbe.SysConstActive(Lx=5, Ly=5)
    state = Swalbe.Sys(sys, "CPU")
    state_act = Swalbe.Sys(sys_act, "CPU")
    state.height .= 1.0
    state_act.height .= 1.0
    # One dim model   
    sys1D = Swalbe.SysConst_1D(L=30)
    sys1D_act = Swalbe.SysConstActive_1D(L=30)
    state1D = Swalbe.Sys(sys1D)
    state1D_act = Swalbe.Sys(sys1D_act)
    state1D.height .= 1.0
    state1D_act.height .= 1.0
    @testset "Slippage" begin
        # No velocities
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), zeros(5,5), 1.0, 1/6)
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state_act, sys_act)
        @test all(fx .== 0.0)
        @test all(fy .== 0.0)
        @test all(state.slipx .== 0.0)
        @test all(state.slipy .== 0.0)
        @test all(state_act.slipx .== 0.0)
        @test all(state_act.slipy .== 0.0)

        # Velocity in x
        Swalbe.slippage!(fx, fy, ones(5,5), fill(0.1,5,5), zeros(5,5), 1.0, 1/6)
        state.velx .= 0.1
        state_act.velx .= 0.1
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state_act, sys_act)
        @test all(fx .== 0.1/11)
        @test all(fy .== 0.0)
        # At some point I have to find out how to use all() with atol
        #@test state.slipx[1,1] ≈ 0.1/11 # atol=1e-8
        @test all(isapprox.(state.slipx, 0.1/11; atol=1e-8))
        @test all(state.slipy .== 0.0)
        #@test state_act.slipx[1,1] ≈ 0.1/11 # atol=1e-8
        @test all(isapprox.(state_act.slipx, 0.1/11; atol=1e-8))
        @test all(state_act.slipy .== 0.0)
        # Velocity in y
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), fill(0.1,5,5), 1.0, 1/6)
        state.velx .= 0.0
        state.vely .= 0.1
        state_act.velx .= 0.0
        state_act.vely .= 0.1
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state_act, sys_act)
        @test all(fx .== 0.0)
        @test all(fy .== 0.1/11)
        @test all(state.slipx .== 0.0)
        #@test state.slipy[1,1] ≈ 0.1/11 # atol=1e-8
        @test all(isapprox.(state.slipy, 0.1/11; atol=1e-8))
        @test all(state_act.slipx .== 0.0)
        #@test state_act.slipy[1,1] ≈ 0.1/11 # atol=1e-8
        @test all(isapprox.(state_act.slipy, 0.1/11, atol=1e-8))
        # Velocity
        Swalbe.slippage!(fx, fy, ones(5,5), fill(-0.1,5,5), fill(0.1,5,5), 1.0, 1/6)
        state.velx .= -0.1
        state.vely .= 0.1
        state_act.velx .= -0.1
        state_act.vely .= 0.1
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state_act, sys_act)
        @test all(fx .== -0.1/11)
        @test all(fy .== 0.1/11)
        # @test state.slipx[1,1] ≈ -0.1/11 
        # @test state.slipy[1,1] ≈ 0.1/11
        # @test state_act.slipx[1,1] ≈ -0.1/11 
        # @test state_act.slipy[1,1] ≈ 0.1/11
        @test all(isapprox.(state.slipx, -0.1/11))
        @test all(isapprox.(state.slipy, 0.1/11))
        @test all(isapprox.(state_act.slipx, -0.1/11))
        @test all(isapprox.(state_act.slipy, 0.1/11))
        # No slip
        Swalbe.slippage!(fx, fy, ones(5,5), fill(-0.1,5,5), fill(0.1,5,5), 0.0, 1/6)
        state.velx .= -0.1
        state.vely .= 0.1
        state_act.velx .= -0.1
        state_act.vely .= 0.1
        sys2 = Swalbe.SysConst(Lx=5,Ly=5,δ=0.0)
        sys2_act = Swalbe.SysConstActive(Lx=5,Ly=5,δ=0.0)
        Swalbe.slippage!(state, sys2)
        Swalbe.slippage!(state_act, sys2_act)
        @test all(fx .== -0.1/2)
        @test all(fy .== 0.1/2)
        # @test state.slipx[1,1] ≈ -0.1/2
        # @test state.slipy[1,1] ≈ 0.1/2
        # @test state_act.slipx[1,1] ≈ -0.1/2
        # @test state_act.slipy[1,1] ≈ 0.1/2
        @test all(isapprox.(state.slipx, -0.1/2))
        @test all(isapprox.(state.slipy, 0.1/2))
        @test all(isapprox.(state_act.slipx, -0.1/2))
        @test all(isapprox.(state_act.slipy, 0.1/2))
    end


    @testset "Slippage 1D" begin
        # No velocities
        Swalbe.slippage!(f1, ones(30), zeros(30), 1.0, 1/6)
        Swalbe.slippage!(state1D, sys1D)
        Swalbe.slippage!(state1D_act, sys1D_act)
        @test all(f1 .== 0.0)
        @test all(state1D.slip .== 0.0)
        @test all(state1D_act.slip .== 0.0)
        # Velocity in x
        Swalbe.slippage!(f1, ones(30), fill(0.1,30), 1.0, 1/6)
        state1D.vel .= 0.1
        state1D.slip .= 0.0
        state1D_act.vel .= 0.1
        state1D_act.slip .= 0.0
        Swalbe.slippage!(state1D, sys1D)
        Swalbe.slippage!(state1D_act, sys1D_act)
        @test all(f1 .== 0.1/11)
        # @test state1D.slip[1] ≈ 0.1/11
        # @test state1D_act.slip[1] ≈ 0.1/11
        @test all(isapprox.(state1D.slip, 0.1/11))
        @test all(isapprox.(state1D_act.slip, 0.1/11))
        # No slip
        Swalbe.slippage!(f1, ones(30), fill(-0.1,30), 0.0, 1/6)
        state1D.vel .= -0.1
        state1D.slip .= 0.0
        state1D_act.vel .= -0.1
        state1D_act.slip .= 0.0
        sys1D2 = Swalbe.SysConst_1D(L=30, δ=0.0)
        sys1D2_act = Swalbe.SysConstActive_1D(L=30, δ=0.0)
        Swalbe.slippage!(state1D, sys1D2)
        Swalbe.slippage!(state1D_act, sys1D2_act)
        @test all(f1 .== -0.1/2)
        # @test state1D.slip[1] ≈ -0.1/2
        # @test state1D_act.slip[1] ≈ -0.1/2
        @test all(isapprox.(state1D.slip, -0.1/2))
        @test all(isapprox.(state1D_act.slip, -0.1/2))
    end

    @testset "Pressure gradient" begin
        sys = Swalbe.SysConst(Lx=5, Ly=5)
        state = Swalbe.Sys(sys, "CPU")    
        state.pressure .= reshape(collect(1.0:25),5,5)
        state.height .= 1.0
        Swalbe.h∇p!(state)

        sys_act = Swalbe.SysConstActive(Lx=5, Ly=5)
        state_act = Swalbe.Sys(sys_act, "CPU")    
        state_act.pressure .= reshape(collect(1.0:25),5,5)
        state_act.height .= 1.0
        Swalbe.h∇p!(state_act)
        
        solx = [-1.5 -1.5 -1.5 -1.5 -1.5;
                 1.0 1.0 1.0 1.0 1.0;
                 1.0 1.0 1.0 1.0 1.0;
                 1.0 1.0 1.0 1.0 1.0;
                -1.5 -1.5 -1.5 -1.5 -1.5]
        soly = [-7.5 5.0 5.0 5.0 -7.5;
                -7.5 5.0 5.0 5.0 -7.5;
                -7.5 5.0 5.0 5.0 -7.5;
                -7.5 5.0 5.0 5.0 -7.5;
                -7.5 5.0 5.0 5.0 -7.5]
        
        @test all(state.h∇px .== solx)
        @test all(state.h∇py .== soly)
        @test all(state_act.h∇px .== solx)
        @test all(state_act.h∇py .== soly)
    end

    @testset "Thermal" begin
        f1 = ones(50,50)
        f2 = ones(50,50)
        vartest = 2*0.01/11
        Swalbe.thermal!(f1, f2, ones(50,50), 0.01, 1/6, 1.0)
        @test mean(f1) ≈ 0.0 atol=1e-2
        @test mean(f2) ≈ 0.0 atol=1e-2
        @test var(f1) ≈ vartest atol=vartest/10
        @test var(f2) ≈ vartest atol=vartest/10
        vartest = 0.2/11
        Swalbe.thermal!(f1, f2, ones(50,50), 0.1, 1/6, 1.0)
        @test mean(f1) ≈ 0.0 atol=1e-2
        @test mean(f2) ≈ 0.0 atol=1e-2
        @test var(f1) ≈ vartest atol=vartest/10
        @test var(f2) ≈ vartest atol=vartest/10
    end
    @testset "Thermal 1D" begin
        f1D = ones(100000)
        vartest = 2*0.01/11
        Swalbe.thermal!(f1D, ones(100000), 0.01, 1/6, 1.0)
        @test mean(f1D) ≈ 0.0 atol=1e-2
        @test var(f1D) ≈ vartest atol=vartest/10
        vartest = 0.2/11
        Swalbe.thermal!(f1D, ones(100000), 0.1, 1/6, 1.0)
        @test mean(f1D) ≈ 0.0 atol=1e-2
        @test var(f1D) ≈ vartest atol=vartest/10
    end

    

    @testset "view four" begin
        dummy = reshape(collect(1:20),5,4)
        d1, d2, d3, d4 = Swalbe.view_four(dummy)
        @test all(d1 .== dummy[:,1])
        @test all(d2 .== dummy[:,2])
        @test all(d3 .== dummy[:,3])
        @test all(d4 .== dummy[:,4])
    end
end