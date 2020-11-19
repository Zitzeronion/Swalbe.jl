@testset "Capillary pressure" begin
    f = reshape(collect(1.0:25),5,5)
    f_float = reshape(collect(1.0f0:25.0f0),5,5)
    res = zeros(5,5)
    @testset "No contact angle" begin
        Swalbe.filmpressure!(res, f, 1.0, 0.0, 3, 2, 0.1, 0.1)
        sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0;
               ]
        for i in eachindex(sol)
            @test res[i] .≈ sol[i] atol=1e-10
        end
    end
    @testset "No height gradient" begin
        nograd = ones(5,5)
        Swalbe.filmpressure!(res, nograd, 1.0, 1/2, 3, 2, 0.1, 0.0)
        for i in eachindex(res)
            @test res[i] .≈ -2(0.1^2-0.1) atol=1e-10
        end
    end
    @testset "Gradient and contact angle" begin
        Swalbe.filmpressure!(res, f, 1.0, 1/2, 3, 2, 0.1, 0.0)
        sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0;
               ]
        println()
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            @test res[i] .≈ sol[i] - 20((0.1/f[i])^3-(0.1/f[i])^2) atol=1e-10
        end
    end
    @testset "Gradient and contact angle Float32" begin
        Swalbe.filmpressure!(res, f_float, 1.0f0, 0.5f0, 3, 2, 0.1f0, 0.0f0)
        sol = [-30.0f0 -5.0f0 -5.0f0 -5.0f0 20f0;
               -25.0f0 0.0f0 0.0f0 0.0f0 25.0f0;
               -25.0f0 0.0f0 0.0f0 0.0f0 25.0f0;
               -25.0f0 0.0f0 0.0f0 0.0f0 25.0f0;
               -20.0f0 5.0f0 5.0f0 5.0f0 30.0f0;
               ]
 
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            @test res[i] .≈ sol[i] - 20((0.1f0/f_float[i])^3-(0.1f0/f_float[i])^2) atol=1e-6
        end
    end
end

@testset "Power function" begin
    argsimple = 5.0
    argsimplefp32 = 5.0f0
    argvec = [2.0 3.0 4.0]
    n1 = 3
    n2 = 6
    
    res0 = Swalbe.power_broad(argsimple, n1)
    @test res0 == 125.0
    res0 = Swalbe.power_broad(argsimple, n2)
    @test res0 = 15625.0
end