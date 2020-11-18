@testset "Finite differences" begin
    f = reshape(collect(1.0:25),5,5)
    # Helps to validate the testing, we use it for DiffEqOperators
    fpad = padarray(f, Pad(:circular,1,1))
    @testset "Gradients" begin
        outputx = zeros(5,5)
        outputy = zeros(5,5)
        @testset "Simple" begin    
            Swalbe.∇f!(outputx, outputy, f)
            @test isa(outputx, Array)
            @test isa(outputy, Array)
            # Analytical results
            solx = zeros(5,5)
            soly = zeros(5,5)
            Dx = CenteredDifference{1}(1,2,1.0,5)
            Dy = CenteredDifference{2}(1,2,1.0,5)
            mul!(solx, Dx, fpad)
            mul!(soly, Dy, fpad)
            # Test them
            @test all(outputx .== solx)
            @test all(outputy .== soly)
        end
        @testset "Four Arguments" begin
            a = fill(0.1, 5,5)
            Swalbe.∇f!(outputx, outputy, f, a)
            @test isa(outputx, Array)
            @test isa(outputy, Array)
            # Analytical results
            solx = zeros(5,5)
            soly = zeros(5,5)
            Dx = CenteredDifference{1}(1,2,1.0,5)
            Dy = CenteredDifference{2}(1,2,1.0,5)
            mul!(solx, Dx, fpad)
            mul!(soly, Dy, fpad)
            # Test them
            @test all(outputx .== 0.1 .* solx)
            @test all(outputy .== 0.1 .* soly)
        end
    end

    @testset "Laplacian" begin
        output = zeros(5,5)
        γ = -1.0
        Swalbe.∇²f!(output, f, γ)
        @test isa(output, Array)
        # println(output)
        sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0;
               ]

        for i in eachindex(sol)
            @test output[i] .≈ sol[i] atol=1e-10
        end
    end
end

