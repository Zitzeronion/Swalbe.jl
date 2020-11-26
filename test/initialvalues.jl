@testset "Initial conditions" begin
    @testset "Fluid Sine Wave" begin
        h = ones(10,10)
        Swalbe.sinewave1D!(h, 2.0, 1, 1, 1)
        @test maximum(h) ≈ 2*(1 + 1) atol= 0.5
        @test minimum(h) ≈ 2*(1 - 1) atol= 0.5
    end
end