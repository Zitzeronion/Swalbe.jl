@testset "Forcings" begin
    fx = zeros(5,5)
    fy = zeros(5,5)
    @testset "Slippage" begin
        # No velocities
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), zeros(5,5), 1.0, 1/6)
        @test all(fx .== 0.0)
        @test all(fy .== 0.0)
        # Velocity in x
        Swalbe.slippage!(fx, fy, ones(5,5), fill(0.1,5,5), zeros(5,5), 1.0, 1/6)
        @test all(fx .== 0.1/11)
        @test all(fy .== 0.0)
        # Velocity in y
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), fill(0.1,5,5), 1.0, 1/6)
        @test all(fx .== 0.0)
        @test all(fy .== 0.1/11)
        # Velocity
        Swalbe.slippage!(fx, fy, ones(5,5), fill(-0.1,5,5), fill(0.1,5,5), 1.0, 1/6)
        @test all(fx .== -0.1/11)
        @test all(fy .== 0.1/11)
        # No slip
        Swalbe.slippage!(fx, fy, ones(5,5), fill(-0.1,5,5), fill(0.1,5,5), 0.0, 1/6)
        @test all(fx .== -0.1/2)
        @test all(fy .== 0.1/2)
    end
end