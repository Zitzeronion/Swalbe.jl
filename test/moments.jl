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
            @test all(i[2] .== 1.0)
        end
    end

    @testset "Velocity" begin
        for i in [state.height, state2.basestate.height, state1D.height, state1D2.basestate.height]
            i .= 0.0
        end
        for i in [f, state.fout, state2.basestate.fout]
            i[:,:,2] .= 0.1
        end
        for i in [f1D, state1D.fout, state1D2.basestate.fout]
            i[:,2] .= 0.1
        end
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(hei, vel, f1D)
        Swalbe.moments!(state)
        Swalbe.moments!(state2)
        Swalbe.moments!(state1D)
        Swalbe.moments!(state1D2)
        @test all(height .== 1.1)
        @test all(hei .== 1.1)
        @test all(velx .== 0.1/1.1)
        @test all(vely .== 0)
        @test all(vel .== 0.1/1.1)

        for i in [state, state2.basestate]
            @test all(i.height .== 1.1)
            @test all(i.velx .== 0.1/1.1)
            @test all(i.vely .== 0.0)
        end
        for i in [state1D, state1D2.basestate]
            @test all(i.height .== 1.1)
            @test all(i.vel .== 0.1/1.1)
        end
        # with f3
        for i in [state.height, state2.basestate.height, state1D.height, state1D2.basestate.height]
            i .= 0.0
        end
        for i in [f, state.fout, state2.basestate.fout]
            i[:,:,3] .= 0.2
        end
        for i in [f1D, state1D.fout, state1D2.basestate.fout]
            i[:,3] .= 0.2
        end
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(hei, vel, f1D)
        Swalbe.moments!(state)
        Swalbe.moments!(state2)
        Swalbe.moments!(state1D)
        Swalbe.moments!(state1D2)
        for i in [state, state2.basestate]
            @test all(i.height .== 1.3)
            @test all(i.velx .== 0.1/1.3)
            @test all(i.vely .== 0.2/1.3)
        end
        for i in [state1D, state1D2.basestate]
            @test all(isapprox.(i.height, 1.3, atol=1e-10))
            @test all(isapprox.(i.vel, -0.1/1.3, atol=1e-10))
        end
        # 2D
        for i in [state.height, state2.basestate.height, state1D.height, state1D2.basestate.height]
            i .= 0.0
        end
        for i in [f, state.fout, state2.basestate.fout]
            i[:,:,4] .= -0.2
        end
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        Swalbe.moments!(state2)
        for i in [state, state2.basestate]
            @test all(i.height .== 1.1)
            @test all(isapprox.(i.velx, 0.3/1.1, atol=1e-10))
            @test all(isapprox.(i.vely, 0.2/1.1, atol=1e-10))
        end
        for i in [state.height, state2.basestate.height, state1D.height, state1D2.basestate.height]
            i .= 0.0
        end
        for i in [f, state.fout, state2.basestate.fout]
            i[:,:,5] .= -0.1
        end
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        Swalbe.moments!(state2)
        for i in [state, state2.basestate]
            @test all(i.height .== 1.0)
            @test all(isapprox.(i.velx, 0.3/1.0, atol=1e-10))
            @test all(isapprox.(i.vely, 0.3/1.0, atol=1e-10))
        end
        for i in [state.height, state2.basestate.height, state1D.height, state1D2.basestate.height]
            i .= 0.0
        end
        for i in [f, state.fout, state2.basestate.fout]
            i[:,:,6] .= 0.1
        end
        Swalbe.moments!(height, velx, vely, f)
        Swalbe.moments!(state)
        Swalbe.moments!(state2)
        for i in [state, state2.basestate]
            @test all(i.height .== 1.1)
            @test all(isapprox.(i.velx, 0.4/1.1, atol=1e-10))
            @test all(isapprox.(i.vely, 0.4/1.1, atol=1e-10))
        end
    end 
    
end