@testset "Moments" begin
    f = zeros(5,5,9)
    f1D = zeros(30,3)
    height = zeros(5,5)
    hei = zeros(30)
    velx = zeros(5,5)
    vely = zeros(5,5)
    vel = zeros(30)
    @testset "No velocity" begin
        f[:,:,1] .= 1.0
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.0)
        f1D[:,1] .= 1.0
        Swalbe.moments!(hei, vel, f1D)
        @test all(hei .== 1.0)
    end
    @testset "Velocity" begin
        # 2D
        f[:,:,1] .= 1.0
        f[:,:,2] .= 0.1
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.1)
        @test all(velx .== 0.1/1.1)
        # 1D
        f1D[:,1] .= 1.0
        f1D[:,2] .= 0.1
        Swalbe.moments!(hei, vel, f1D)
        @test all(hei .== 1.1)
        @test all(vel .== 0.1/1.1)
        # 2D
        f[:,:,3] .= 0.2
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.3)
        @test all(velx .== 0.1/1.3)
        @test all(vely .== 0.2/1.3)
        # 1D
        f1D[:,3] .= 0.2
        Swalbe.moments!(hei, vel, f1D)
        @test all(hei .== 1.3)
        @test all(vel .== -0.1/1.3)
        # 2D
        f[:,:,4] .= -0.2
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.1)
        @test all(velx .≈ 0.3/1.1)
        @test all(vely .== 0.2/1.1)
        f[:,:,5] .= -0.1
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.0)
        @test all(velx .≈ 0.3/1.0)
        @test all(vely .≈ 0.3/1.0)
        f[:,:,6] .= 0.1
        Swalbe.moments!(height, velx, vely, f)
        @test all(height .== 1.1)
        @test all(velx .≈ 0.4/1.1)
        @test all(vely .≈ 0.4/1.1)
    end 
end