@testset "Initial conditions" begin
    @testset "Fluid Sine Wave" begin
        h = ones(10,10)
        Swalbe.sinewave1D!(h, 2.0, 1, 1, 1)
        @test maximum(h) ≈ 2*(1 + 1) atol= 0.5
        @test minimum(h) ≈ 2*(1 - 1) atol= 0.5
    end

    @testset "Random interface" begin
        h = ones(10,10)
        Swalbe.randinterface!(h, 2.0, 0.1)
        @test maximum(h) <= 2.2
        @test minimum(h) >= 1.8
    end

    @testset "Single droplet" begin
        rad = 50
        θ = 1/3
        center = (50,50)
        height = Swalbe.singledroplet(ones(100,100), rad, θ, center)
        @test isa(height, Array{Float64,2})
        @test size(height) == (100, 100)
        @test findmax(height)[1] == rad * (1 - cospi(θ))
        @test findmax(height)[2] == CartesianIndex(center[1],center[2])
    end
end
@testset "Contact angle patterns" begin
    @testset "Default elliptical patch" begin
        theta, P = Swalbe.ellipsepattern(ones(100,100), 1/9) # with kw: center=(Int(lx/2), Int(ly/2)),δₐ=1/36, a=10, b=5
        @test size(theta) == (100, 100)
        @test theta[50,50] == 1/9 + 1/36
        @test theta[1,1] == 1/9
    end
    @testset "Off center elliptical patch" begin
        theta, P = Swalbe.ellipsepattern(ones(100,100), 1/9, center=(20,20)) # with kw: center=(Int(lx/2), Int(ly/2)),δₐ=1/36, a=10, b=5
        @test size(theta) == (100, 100)
        @test theta[20,20] == 1/9 + 1/36
        @test theta[100,100] == 1/9
    end
    @testset "Negative contrast elliptical patch" begin
        theta, P = Swalbe.ellipsepattern(ones(100,100), 1/9, δₐ=-1/36) # with kw: center=(Int(lx/2), Int(ly/2)),δₐ=1/36, a=10, b=5
        @test size(theta) == (100, 100)
        @test theta[50,50] == 1/9 - 1/36
        @test theta[100,100] == 1/9
    end
    @testset "Larger elliptical patch" begin
        theta, P = Swalbe.ellipsepattern(ones(100,100), 1/9, a=10, b=20) # with kw: center=(Int(lx/2), Int(ly/2)),δₐ=1/36, a=10, b=5
        @test size(theta) == (100, 100)
        @test theta[40,50] == 1/9 + 1/36
        @test theta[70,50] == 1/9 + 1/36
        @test theta[100,100] == 1/9
    end
    @testset "Polygon patch" begin
        lx, ly = 100, 100  
        @testset "Equilateral Triangle" begin
            θₙ, P = Swalbe.trianglepattern(ones(100, 100), 1/9)
            @test all(θₙ .> 0.0)
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(lx/2), Int(ly/2)] == 1/9 + 1/36 
            θₙ, P = Swalbe.trianglepattern(ones(100, 100), 1/9, center=(lx/3, 2*ly/3))
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(round(lx/3)), Int(round(2*ly/3))] == 1/9 # Because the side became to large
            @test θₙ[30, 20] == 1/9 + 1/36 # If too large stick to origin!
            θₙ, P = Swalbe.trianglepattern(ones(100, 100), 1/9, δₐ = -1/36)
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(lx/2), Int(ly/2)] == 1/9 - 1/36 
            θₙ, P = Swalbe.trianglepattern(ones(100, 100), 1/9, side = 30) # Smaller triangle
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(lx/2), Int(ly/2)] == 1/9 + 1/36
            @test θₙ[50, 60] == 1/9 + 1/36 
        end

        @testset "Box pattern" begin
            θₙ, P = Swalbe.boxpattern(ones(100, 100), 1/9)
            @test all(θₙ .> 0.0)
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(lx/2), Int(ly/2)] == 1/9 + 1/36 
            θₙ, P = Swalbe.boxpattern(ones(100, 100), 1/9, center=(lx/3, 2*ly/3))
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(round(2*ly/3)), Int(round(lx/3))] == 1/9 + 1/36
            θₙ, P = Swalbe.boxpattern(ones(100, 100), 1/9, δₐ = -1/36)
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(lx/2), Int(ly/2)] == 1/9 - 1/36 
            θₙ, P = Swalbe.boxpattern(ones(100, 100), 1/9, side = 30) # Smaller triangle
            @test θₙ[1,1] == 1/9 # the default value
            @test θₙ[Int(lx/2), Int(ly/2)] == 1/9 + 1/36
            @test θₙ[Int(lx/2) + Int(30/2), Int(ly/2) + Int(30/2)] == 1/9 + 1/36 
        end
    end
end