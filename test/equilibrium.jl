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