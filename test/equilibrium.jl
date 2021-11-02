@testset "Equilibria" begin
    @testset "Nothing" begin
        feq = zeros(5,5,9)
        Swalbe.equilibrium!(feq, zeros(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.0)
        @test all(feq[:,:,:] .== 0.0)
        Swalbe.equilibrium!(feq, zeros(5,5), zeros(5,5), zeros(5,5), zeros(5,5))
        @test all(feq[:,:,:] .== 0.0)
    end
    @testset "Density only" begin
        feq = zeros(5,5,9)
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.0)
        @test all(feq[:,:,1] .== 1.0)
        @test all(feq[:,:,2:9] .== 0.0)
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5))
        @test all(feq[:,:,1] .== 1.0)
        @test all(feq[:,:,2:9] .== 0.0)
    end
    @testset "Density and gravity" begin
        feq = zeros(5,5,9)
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.1)
        @test all(feq[:,:,1] .== 1.0 - 1/12)
        @test all(feq[:,:,2:5] .≈ 1/9 * 1.5 * 0.1)
        @test all(feq[:,:,6:9] .≈ 1/36 * 1.5 * 0.1)
    end
    @testset "Density and velocity" begin
        feq = zeros(5,5,9)
        # with gravity
        Swalbe.equilibrium!(feq, ones(5,5), fill(0.1,5,5), fill(-0.1,5,5), zeros(5,5), 0.0)
        @test all(feq[:,:,1] .== 1.0 - 2/3 * 0.02)
        @test all(feq[:,:,2] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,3] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,4] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,5] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,6] .≈  1/36 * (-3/2 * 0.02))
        @test all(feq[:,:,7] .≈  1/36 * (3 * -0.2  + 4.5 * 0.2^2 - 3/2 * 0.02))
        @test all(feq[:,:,8] .≈  1/36 * (-3/2 * 0.02))
        @test all(feq[:,:,9] .≈  1/36 * (3 * 0.2 + 4.5 * 0.2^2 - 3/2 * 0.02))
        # without gravity
        Swalbe.equilibrium!(feq, ones(5,5), fill(0.1,5,5), fill(-0.1,5,5), zeros(5,5))
        @test all(feq[:,:,1] .== 1.0 - 2/3 * 0.02)
        @test all(feq[:,:,2] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,3] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,4] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,5] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,6] .≈  1/36 * (-3/2 * 0.02))
        @test all(feq[:,:,7] .≈  1/36 * (3 * -0.2  + 4.5 * 0.2^2 - 3/2 * 0.02))
        @test all(feq[:,:,8] .≈  1/36 * (-3/2 * 0.02))
        @test all(feq[:,:,9] .≈  1/36 * (3 * 0.2 + 4.5 * 0.2^2 - 3/2 * 0.02))
        # using the state struct
        state = Swalbe.Sys(Swalbe.SysConst(Lx=5, Ly=5), "CPU")
        state.feq .= 0.0
        state.height .= 1.0 
        state.velx .= 0.1 
        state.vely .= -0.1
        Swalbe.equilibrium!(state)
        @test all(state.feq[:,:,1] .== 1.0 - 2/3 * 0.02)
        @test all(state.feq[:,:,2] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(state.feq[:,:,3] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(state.feq[:,:,4] .≈  1/9 * (3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(state.feq[:,:,5] .≈  1/9 * (3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(state.feq[:,:,6] .≈  1/36 * (-3/2 * 0.02))
        @test all(state.feq[:,:,7] .≈  1/36 * (3 * -0.2  + 4.5 * 0.2^2 - 3/2 * 0.02))
        @test all(state.feq[:,:,8] .≈  1/36 * (-3/2 * 0.02))
        @test all(state.feq[:,:,9] .≈  1/36 * (3 * 0.2 + 4.5 * 0.2^2 - 3/2 * 0.02))

    end
    @testset "Density and gravity and velocity" begin
        feq = zeros(5,5,9)
        Swalbe.equilibrium!(feq, ones(5,5), fill(0.1,5,5), fill(-0.1,5,5), zeros(5,5), 0.1)
        @test all(feq[:,:,1] .== 1.0 - 1/12 - 2/3 * 0.02)
        @test all(feq[:,:,2] .≈  1/9 * (1.5 * 0.1 + 3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,3] .≈  1/9 * (1.5 * 0.1 + 3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,4] .≈  1/9 * (1.5 * 0.1 + 3 * -0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,5] .≈  1/9 * (1.5 * 0.1 + 3 * 0.1 + 4.5 * 0.01 - 3/2 * 0.02))
        @test all(feq[:,:,6] .≈  1/36 * (1.5 * 0.1 + -3/2 * 0.02))
        @test all(feq[:,:,7] .≈  1/36 * (1.5 * 0.1 + 3 * -0.2  + 4.5 * 0.2^2 - 3/2 * 0.02))
        @test all(feq[:,:,8] .≈  1/36 * (1.5 * 0.1 + -3/2 * 0.02))
        @test all(feq[:,:,9] .≈  1/36 * (1.5 * 0.1 + 3 * 0.2 + 4.5 * 0.2^2 - 3/2 * 0.02))
    end
end

@testset "Equilibria 1D" begin
    @testset "Nothing" begin
        feq = zeros(30,3)
        Swalbe.equilibrium!(feq, zeros(30), zeros(30), 0.0)
        @test all(feq[:,:] .== 0.0)
        Swalbe.equilibrium!(feq, zeros(30), zeros(30))
        @test all(feq[:,:] .== 0.0)
    end
    @testset "Density only" begin
        feq = zeros(30,3)
        Swalbe.equilibrium!(feq, ones(30), zeros(30), 0.0)
        @test all(feq[:,1] .== 1.0)
        @test all(feq[:,2:3] .== 0.0)
        Swalbe.equilibrium!(feq, ones(30), zeros(30))
        @test all(feq[:,1] .== 1.0)
        @test all(feq[:,2:3] .== 0.0)
    end
    @testset "Density and gravity" begin
        feq = zeros(30,3)
        Swalbe.equilibrium!(feq, ones(30), zeros(30), 0.1)
        @test all(feq[:,1] .== 1.0 - 0.05)
        @test all(feq[:,2:3] .== 0.025)
        
    end
    @testset "Density and velocity" begin
        feq = zeros(30,3)
        Swalbe.equilibrium!(feq, ones(30), fill(0.1,30), 0.0)
        @test all(feq[:,1] .== 1.0 - 0.01)
        @test all(feq[:,2] .≈  0.5 * 0.1 + 0.5 * 0.01)
        @test all(feq[:,3] .≈  -0.5 * 0.1 + 0.5 * 0.01)
        
        Swalbe.equilibrium!(feq, ones(30), fill(0.1,30))
        @test all(feq[:,1] .== 1.0 - 0.01)
        @test all(feq[:,2] .≈  0.5 * 0.1 + 0.5 * 0.01)
        @test all(feq[:,3] .≈  -0.5 * 0.1 + 0.5 * 0.01)

        state = Swalbe.Sys(Swalbe.SysConst_1D(L=30))
        state.feq .= 0.0
        state.height .= 1.0
        state.vel .= 0.1
        Swalbe.equilibrium!(state)
        @test all(state.feq[:,1] .== 1.0 - 0.01)
        @test all(state.feq[:,2] .≈  0.5 * 0.1 + 0.5 * 0.01)
        @test all(state.feq[:,3] .≈  -0.5 * 0.1 + 0.5 * 0.01)
        
    end
    @testset "Density and gravity and velocity" begin
        feq = zeros(30,3)
        Swalbe.equilibrium!(feq, ones(30), fill(0.1,30), 0.1)
        @test all(feq[:,1] .== 1.0 - 0.5 * 0.1 - 0.01)
        @test all(feq[:,2] .≈  0.025 + 0.5 * 0.1 + 0.5 * 0.01)
        @test all(feq[:,3] .≈  0.025 - 0.5 * 0.1 + 0.5 * 0.01)
    end
end