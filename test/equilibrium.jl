@testset "Equilibria" begin
    @testset "Nothing" begin
    feq = zeros(5,5,9)
    Swalbe.equilibrium!(feq, zeros(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.0)
    @test all(feq[:,:,:] .== 0.0)
    end
    @testset "Just density" begin
        feq = zeros(5,5,9)
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.0)
        @test all(feq[:,:,1] .== 1.0)
        @test all(feq[:,:,2:9] .== 0.0)
    end
    @testset "Density and gravity" begin
        feq = zeros(5,5,9)
        Swalbe.equilibrium!(feq, ones(5,5), zeros(5,5), zeros(5,5), zeros(5,5), 0.1)
        @test all(feq[:,:,1] .== 1.0 - 1/12)
        @test all(feq[:,:,2:5] .≈ 1/9 * 1.5 * 0.1)
        println("Result $(feq[1,1,2]) and $(1/9 * 1.5 * 0.1)")
        @test all(feq[:,:,6:9] .≈ 1/36 * 1.5 * 0.1)
    end
end