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