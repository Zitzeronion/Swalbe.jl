@testset "Forcings" begin
    fx = zeros(5,5)
    fy = zeros(5,5)
    @testset "Slippage" begin
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), zeros(5,5), 1.0, 1/6)
        @test all(fx .== 0.0)
        @test all(fy .== 0.0)
    end
end