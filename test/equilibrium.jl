@testset "Equilibria" begin
    feq = zeros(5,5,9)
    sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs())
    state = Swalbe.Sys(sys, "CPU")
    state2 = Swalbe.Sys(sys, "CPU", kind="thermal")
    @testset "Nothing" begin
        state.height .= 0.0
        state2.basestate.height .= 0.0
        Swalbe.equilibrium!(state, sys, g=0.0)
        Swalbe.equilibrium!(state2, sys, g=0.0)
        Swalbe.equilibrium!(feq, zeros(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.0)
        for i in [feq, state.feq, state2.basestate.feq]
            @test all(i[:,:,:] .== 0.0)
        end
    end
    
    @testset "Density only" begin
        state.height .= 1.0
        state2.basestate.height .= 1.0
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.0)
        Swalbe.equilibrium!(state, sys, g=0.0)
        Swalbe.equilibrium!(state2, sys, g=0.0)
        for i in [feq, state.feq, state2.basestate.feq]
            @test all(i[:,:,1] .== 1.0)
            @test all(i[:,:,2:9] .== 0.0)
        end
    end
    
    @testset "Density and gravity" begin
        state.height .= 1.0
        state2.basestate.height .= 1.0
        Swalbe.equilibrium!(state, sys, g=0.1)
        Swalbe.equilibrium!(state2, sys, g=0.1)
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.1)
        for i in [feq, state.feq, state2.basestate.feq]
            @test all(i[:,:,1] .== 1.0 - 1/12)
            @test all(i[:,:,2:5] .≈ 1/9 * 1.5 * 0.1)
            @test all(i[:,:,6:9] .≈ 1/36 * 1.5 * 0.1)
        end
    end
    
    @testset "Density and velocity" begin
        for i in [state, state2.basestate]
            i.feq .= 0.0
            i.height .= 1.0 
            i.velx .= 0.1 
            i.vely .= -0.1
        end
        Swalbe.equilibrium!(feq, ones(5,5), fill(0.1,5,5), fill(-0.1,5,5), zeros(5,5), 0.0)
        Swalbe.equilibrium!(state, sys)
        Swalbe.equilibrium!(state2, sys)
        for i in [feq, state.feq, state2.basestate.feq]
            @test all(i[:,:,1] .== 1.0 - 2/3 * 0.02)
            @test all(i[:,:,2] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,3] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,4] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,5] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,6] .≈  1/36 * (-3/2 * 0.02))
            @test all(i[:,:,7] .≈  1/36 * (3 * -0.2  + 4.5 * 0.2^2 - 3/2 * 0.02))
            @test all(i[:,:,8] .≈  1/36 * (-3/2 * 0.02))
            @test all(i[:,:,9] .≈  1/36 * (3 * 0.2 + 4.5 * 0.2^2 - 3/2 * 0.02))
        end
    end
    
    @testset "Density and gravity and velocity" begin
        Swalbe.equilibrium!(feq, ones(5,5), fill(0.1,5,5), fill(-0.1,5,5), zeros(5,5), 0.1)
        Swalbe.equilibrium!(state, sys, g=0.1)
        Swalbe.equilibrium!(state2, sys, g=0.1)
        for i in [feq, state.feq, state2.basestate.feq]
            @test all(i[:,:,1] .== 1.0 - 1/12 - 2/3 * 0.02)
            @test all(i[:,:,2] .≈  1/9 * (1.5 * 0.1 + 3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,3] .≈  1/9 * (1.5 * 0.1 + 3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,4] .≈  1/9 * (1.5 * 0.1 + 3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,5] .≈  1/9 * (1.5 * 0.1 + 3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
            @test all(i[:,:,6] .≈  1/36 * (1.5 * 0.1 + -3/2 * 0.02))
            @test all(i[:,:,7] .≈  1/36 * (1.5 * 0.1 + 3 * -0.2  + 4.5 * 0.2^2 - 3/2 * 0.02))
            @test all(i[:,:,8] .≈  1/36 * (1.5 * 0.1 + -3/2 * 0.02))
            @test all(i[:,:,9] .≈  1/36 * (1.5 * 0.1 + 3 * 0.2 + 4.5 * 0.2^2 - 3/2 * 0.02))
        end
    end
    
end

@testset "Equilibria 1D" begin
    feq = zeros(30,3)
    sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs())
    st1 = Swalbe.Sys(sys)
    st2 = Swalbe.Sys(sys, kind="thermal")
    @testset "Nothing" begin
        st1.height .= 0.0
        st2.basestate.height .= 0.0
        Swalbe.equilibrium!(feq, zeros(30), zeros(30), 0.0)
        Swalbe.equilibrium!(st1, sys, g=0.0)
        Swalbe.equilibrium!(st2, sys, g=0.0)
        for i in [feq, st1.feq, st2.basestate.feq]
            @test all(i[:,:] .== 0.0)
        end
    end
    
    @testset "Density only" begin
        st1.height .= 1.0
        st2.basestate.height .= 1.0
        Swalbe.equilibrium!(st1, sys, g=0.0)
        Swalbe.equilibrium!(st2, sys, g=0.0)
        Swalbe.equilibrium!(feq, ones(30), zeros(30), 0.0)
        for i in [feq, st1.feq, st2.basestate.feq]
            @test all(i[:,1] .== 1.0)
            @test all(i[:,2:3] .== 0.0)
        end
    end
    
    @testset "Density and gravity" begin
        Swalbe.equilibrium!(feq, ones(30), zeros(30), 0.1)
        Swalbe.equilibrium!(st1, sys, g=0.1)
        Swalbe.equilibrium!(st2, sys, g=0.1)
        for i in [feq, st1.feq, st2.basestate.feq]
            @test all(i[:,1] .== 1.0 - 0.05)
            @test all(i[:,2:3] .== 0.025)
        end
    end

    @testset "Density and velocity" begin
        for i in [st1, st2.basestate]
            i.vel .= 0.1
        end
        Swalbe.equilibrium!(feq, ones(30), fill(0.1,30), 0.0)
        Swalbe.equilibrium!(st1, sys, g=0.0)
        Swalbe.equilibrium!(st2, sys, g=0.0)
        for i in [feq, st1.feq, st2.basestate.feq]
            @test all(i[:,1] .== 1.0 - 0.01)
            @test all(i[:,2] .≈  0.5 * 0.1 + 0.5 * 0.01)
            @test all(i[:,3] .≈  -0.5 * 0.1 + 0.5 * 0.01)
        end
    end
    
    @testset "Density and gravity and velocity" begin
        Swalbe.equilibrium!(st1, sys, g=0.1)
        Swalbe.equilibrium!(st2, sys, g=0.1)
        Swalbe.equilibrium!(feq, ones(30), fill(0.1,30), 0.1)
        for i in [feq, st1.feq, st2.basestate.feq]
            @test all(feq[:,1] .== 1.0 - 0.5 * 0.1 - 0.01)
            @test all(feq[:,2] .≈  0.025 + 0.5 * 0.1 + 0.5 * 0.01)
            @test all(feq[:,3] .≈  0.025 - 0.5 * 0.1 + 0.5 * 0.01)
        end
    end

end
