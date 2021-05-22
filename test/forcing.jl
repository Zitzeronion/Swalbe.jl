@testset "Forcings" begin
    fx = zeros(5,5)
    fy = zeros(5,5)
    f1 = zeros(30)
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
    @testset "Slippage 1D" begin
        # No velocities
        Swalbe.slippage!(f1, ones(30), zeros(30), 1.0, 1/6)
        @test all(f1 .== 0.0)
        # Velocity in x
        Swalbe.slippage!(f1, ones(30), fill(0.1,30), 1.0, 1/6)
        @test all(f1 .== 0.1/11)
        # No slip
        Swalbe.slippage!(f1, ones(30), fill(-0.1,30), 0.0, 1/6)
        @test all(f1 .== -0.1/2)
    end

    @testset "Thermal" begin
        f1 = ones(50,50)
        f2 = ones(50,50)
        vartest = 2*0.01/11
        Swalbe.thermal!(f1, f2, ones(50,50), 0.01, 1/6, 1.0)
        @test mean(f1) ≈ 0.0 atol=1e-2
        @test mean(f2) ≈ 0.0 atol=1e-2
        @test var(f1) ≈ vartest atol=vartest/10
        @test var(f2) ≈ vartest atol=vartest/10
        vartest = 0.2/11
        Swalbe.thermal!(f1, f2, ones(50,50), 0.1, 1/6, 1.0)
        @test mean(f1) ≈ 0.0 atol=1e-2
        @test mean(f2) ≈ 0.0 atol=1e-2
        @test var(f1) ≈ vartest atol=vartest/10
        @test var(f2) ≈ vartest atol=vartest/10
    end
    @testset "Thermal 1D" begin
        f1D = ones(100000)
        vartest = 2*0.01/11
        Swalbe.thermal!(f1D, ones(100000), 0.01, 1/6, 1.0)
        @test mean(f1D) ≈ 0.0 atol=1e-2
        @test var(f1D) ≈ vartest atol=vartest/10
        vartest = 0.2/11
        Swalbe.thermal!(f1D, ones(100000), 0.1, 1/6, 1.0)
        @test mean(f1D) ≈ 0.0 atol=1e-2
        @test var(f1D) ≈ vartest atol=vartest/10
    end

    @testset "rho update" begin
        @testset "Constant fields" begin
            rho = ones(25)
            height = ones(25)
            output = zeros(25)
            Swalbe.update_rho!(rho, output, height, zeros(25,2), zeros(25,4))
            @test all(rho .== 1)
        end
    end

    @testset "view four" begin
        dummy = reshape(collect(1:20),5,4)
        d1, d2, d3, d4 = Swalbe.view_four(dummy)
        @test all(d1 .== dummy[:,1])
        @test all(d2 .== dummy[:,2])
        @test all(d3 .== dummy[:,3])
        @test all(d4 .== dummy[:,4])
    end
end