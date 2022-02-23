@testset "Moments" begin
    f = zeros(5,5,9)
    f1D = zeros(30,3)
    height = zeros(5,5)
    hei = zeros(30)
    velx = zeros(5,5)
    vely = zeros(5,5)
    vel = zeros(30)
    state = Swalbe.Sys(Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs()), "CPU")
    state1D = Swalbe.Sys(Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs())) 
    state2 = Swalbe.Sys(Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs()), "CPU", kind="thermal")
    state1D2 = Swalbe.Sys(Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs()), kind="thermal") 
    
    @testset "No velocity" begin
        for i in [state.height, state2.basestate.height, state1D.height, state1D2.basestate.height]
            i .= 0.0
        end
        for i in [f, state.fout, state2.basestate.fout]
            i[:,:,1] .= 1.0
        end
        for i in [f1D, state1D.fout, state1D2.basestate.fout]
            i[:,1] .= 1.0
        end
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(hei, vel, f1D)
        Swalbe.moments!(state)
        Swalbe.moments!(state2)
        Swalbe.moments!(state1D)
        Swalbe.moments!(state1D2)
        for i in enumerate([height, state.height, state2.basestate.height, state1D.height, state1D2.basestate.height])
            # println(i[1])
            @test all(i[2] .== 1.0)
        end
    end
    #=
    @testset "Velocity" begin
        # 2D
        f[:,:,1] .= 1.0
        f[:,:,2] .= 0.1
        state.fout[:,:,2] .= 0.1
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        @test all(height .== 1.1)
        @test all(velx .== 0.1/1.1)
        @test all(state.height .== 1.1)
        @test all(state.velx .== 0.1/1.1)
        # 1D
        f1D[:,1] .= 1.0
        f1D[:,2] .= 0.1
        state1D.fout[:,2] .= 0.1
        Swalbe.moments!(hei, vel, f1D)
        @test all(hei .== 1.1)
        @test all(vel .== 0.1/1.1)
        Swalbe.moments!(state1D)
        @test all(state1D.height .== 1.1)
        @test all(state1D.vel .== 0.1/1.1)
        # 2D
        f[:,:,3] .= 0.2
        state.fout[:,:,3] .= 0.2
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        @test all(height .== 1.3)
        @test all(velx .== 0.1/1.3)
        @test all(vely .== 0.2/1.3)
        @test all(state.height .== 1.3)
        @test all(state.velx .== 0.1/1.3)
        @test all(state.vely .== 0.2/1.3)
        # 1D
        f1D[:,3] .= 0.2
        state1D.fout[:,3] .= 0.2
        Swalbe.moments!(hei, vel, f1D)
        Swalbe.moments!(state1D)
        @test all(hei .== 1.3)
        @test all(vel .== -0.1/1.3)
        @test state1D.height[1] ≈ 1.3
        @test state1D.vel[1] ≈ -0.1/1.3
        # 2D
        f[:,:,4] .= -0.2
        state.fout[:,:,4] .= -0.2
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        @test all(height .== 1.1)
        @test all(velx .≈ 0.3/1.1)
        @test all(vely .== 0.2/1.1)
        @test all(state.height .== 1.1)
        @test all(state.velx .≈ 0.3/1.1)
        @test all(state.vely .== 0.2/1.1)
        f[:,:,5] .= -0.1
        state.fout[:,:,5] .= -0.1
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        @test all(height .== 1.0)
        @test all(velx .≈ 0.3/1.0)
        @test all(vely .≈ 0.3/1.0)
        @test all(state.height .== 1.0)
        @test all(state.velx .≈ 0.3/1.0)
        @test all(state.vely .≈ 0.3/1.0)
        f[:,:,6] .= 0.1
        state.fout[:,:,6] .= 0.1
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        @test all(height .== 1.1)
        @test all(velx .≈ 0.4/1.1)
        @test all(vely .≈ 0.4/1.1)
        @test all(state.height .== 1.1)
        @test all(state.velx .≈ 0.4/1.1)
        @test all(state.vely .≈ 0.4/1.1)
    end 
    =#
end