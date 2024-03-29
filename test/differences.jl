@testset "Finite differences" begin
    f = reshape(collect(1.0:25),5,5)
    @testset "Gradients" begin
        outputx = zeros(5,5)
        outputy = zeros(5,5)
        @testset "Simple" begin    
            Swalbe.∇f!(outputx, outputy, f)
            @test isa(outputx, Array)
            @test isa(outputy, Array)
            #= 
            Analytical results, tested with DiffEqOperators
            Dx = CenteredDifference{1}(1,2,1.0,5)
            Dy = CenteredDifference{2}(1,2,1.0,5)
            fpad = padarray(f, Pad(:circular, 1, 1))
            mul!(solx, Dx, fpad)
            mul!(soly, Dy, fpad)
            =#
            solx = [-1.5 -1.5 -1.5 -1.5 -1.5;
                     1.0 1.0 1.0 1.0 1.0;
                     1.0 1.0 1.0 1.0 1.0;
                     1.0 1.0 1.0 1.0 1.0;
                    -1.5 -1.5 -1.5 -1.5 -1.5]

            soly = [-7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5]
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
            solx = [-1.5 -1.5 -1.5 -1.5 -1.5;
                     1.0 1.0 1.0 1.0 1.0;
                     1.0 1.0 1.0 1.0 1.0;
                     1.0 1.0 1.0 1.0 1.0;
                    -1.5 -1.5 -1.5 -1.5 -1.5]

            soly = [-7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5]
            # Test them
            @test all(outputx .== 0.1 .* solx)
            @test all(outputy .== 0.1 .* soly)
        end
        @testset "Five Arguments" begin
            a = fill(0.1, 5,5)
            b = zeros(5,5,8)
            Swalbe.∇f!(outputx, outputy, f, b, a)
            @test isa(outputx, Array)
            @test isa(outputy, Array)
            # Analytical results
            solx = [-1.5 -1.5 -1.5 -1.5 -1.5;
                     1.0 1.0 1.0 1.0 1.0;
                     1.0 1.0 1.0 1.0 1.0;
                     1.0 1.0 1.0 1.0 1.0;
                    -1.5 -1.5 -1.5 -1.5 -1.5]

            soly = [-7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5;
                    -7.5 5.0 5.0 5.0 -7.5]
            # Test them
            @test all(outputx .== 0.1 .* solx)
            @test all(outputy .== 0.1 .* soly)
        end
    end

    f1 = collect(1.0:25)
    @testset "Gradients 1D" begin
        output = zeros(25)
        dummy = zeros(25,2)
        @testset "Simple" begin    
            Swalbe.∇f!(output, f1, dummy, ones(25))
            @test isa(output, Vector)

            sol = ones(25)
            sol[1] = -11.5 
            sol[end] = -11.5 
            # Test them
            @test all(output .== sol)

            Swalbe.∇f!(output, f1, dummy)
            @test isa(output, Vector)

            sol = ones(25)
            sol[1] = -11.5 
            sol[end] = -11.5 
            # Test them
            @test all(output .== sol)
        end
        
    end

    @testset "Laplacian" begin
        output = zeros(5,5)
        γ = -1.0
        Swalbe.∇²f!(output, f, γ)
        @test isa(output, Array)
        #= 
        Analytical results, tested with DiffEqOperators
        Dxx = CenteredDifference{1}(2,2,1.0,5)
        Dyy = CenteredDifference{2}(2,2,1.0,5)
        fpad = padarray(f, Pad(:circular, 1, 1))
        A = Dxx + Dyy
        mul!(sol, A, fpad)
        =#
        sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0];
        
        @test all(isapprox.(output, sol; atol=1e-10))
    end

    @testset "Laplacian 1D" begin
    output = zeros(25)
    dummy = zeros(25,2)
    input = collect(1:25.0)
    Swalbe.∇²f!(output, input, dummy)
    @test isa(output, Vector)
    
    sol = zeros(25)
    sol[1] = 25;
    sol[end] = -25;
    
    @test all(output .== sol)
    
end

    @testset "Viewneighbors" begin
        f = reshape(collect(1.0:5*5*8),5,5,8)
        f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewneighbors(f)
        allviews = [f1, f2, f3, f4, f5, f6, f7, f8]
        for (index, value) in enumerate(allviews)
            @test all(value .== f[:,:,index])
        end
    end

    @testset "Viewneighbors 1D" begin
        f = reshape(collect(1.0:20),10,2)
        f1, f2 = Swalbe.viewneighbors_1D(f)
        allviews = [f1, f2]
        for (index, value) in enumerate(allviews)
            @test all(value .== f[:,index])
        end
    end
end


