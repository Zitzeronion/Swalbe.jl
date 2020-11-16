@testset "Finite differences" begin
    @testset "Gradients" begin
        outputx = zeros(5,5)
        outputy = zeros(5,5)
        f = reshape(collect(1.0:25),5,5)
        Swalbe.∇f!(outputx, outputy, f)
        @test isa(outputx, Array)
        @test isa(outputy, Array)
        println(outputx, outputy)
        # Analytical results
        solx = [1.5 1.5 1.5 1.5 1.5;
               -1.0 -1.0 -1.0 -1.0 -1.0;
               -1.0 -1.0 -1.0 -1.0 -1.0;
               -1.0 -1.0 -1.0 -1.0 -1.0;
                1.5 1.5 1.5 1.5 1.5]
        soly = [7.5 -5.0 -5.0 -5.0 7.5; 
                7.5 -5.0 -5.0 -5.0 7.5; 
                7.5 -5.0 -5.0 -5.0 7.5; 
                7.5 -5.0 -5.0 -5.0 7.5; 
                7.5 -5.0 -5.0 -5.0 7.5]
        # Test them
        @test all(outputx .== solx)
        @test all(outputy .== soly)
    end


end

