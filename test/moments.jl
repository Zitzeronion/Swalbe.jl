@testset "Moments" begin
    f = zeros(5,5,9)
    height = zeros(5,5)
    velx = zeros(5,5)
    vely = zeros(5,5)
    @testset "No velocity" begin
        f[:,:,1] .= 1.0
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.0)
    end
    @testset "Velocity" begin
        f[:,:,1] .= 1.0
        f[:,:,2] .= 0.1
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.1)
        @test all(velx .== 0.1/1.1)
        f[:,:,3] .= 0.2
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.3)
        @test all(velx .== 0.1/1.3)
        @test all(vely .== 0.2/1.3)
        f[:,:,4] .= -0.2
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.1)
        @test all(velx .≈ 0.3/1.1)
        @test all(vely .== 0.2/1.1)
    end 
end