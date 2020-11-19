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
            @test res[i] .≈ 2(0.1^2-0.1) atol=1e-10
        end
    end
    @testset "Gradient and contact angle" begin
        Swalbe.filmpressure!(res, f, 1.0, 1/2, 3, 2, 0.1, 0.0)
        sol = [30.0 5.0 5.0 5.0 -20;
               25.0 0.0 0.0 0.0 -25.0;
               25.0 0.0 0.0 0.0 -25.0;
               25.0 0.0 0.0 0.0 -25.0;
               20.0 -5.0 -5.0 -5.0 -30.0]
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            @test res[i] .≈ -1(sol[i] - 20((0.1/f[i])^3-(0.1/f[i])^2)) atol=1e-10
        end
    end
    @testset "Gradient and contact angle Float32" begin
        Swalbe.filmpressure!(res, f_float, 1.0f0, 0.5f0, 3, 2, 0.1f0, 0.0f0)
        sol = [30.0f0 5.0f0 5.0f0 5.0f0 -20f0;
               25.0f0 0.0f0 0.0f0 0.0f0 -25.0f0;
               25.0f0 0.0f0 0.0f0 0.0f0 -25.0f0;
               25.0f0 0.0f0 0.0f0 0.0f0 -25.0f0;
               20.0f0 -5.0f0 -5.0f0 -5.0f0 -30.0f0]
 
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            @test res[i] .≈ -1(sol[i] - 20((0.1f0/f_float[i])^3-(0.1f0/f_float[i])^2)) atol=1e-6
        end
    end
end

@testset "Power function" begin
    argsimple = 5.0
    argsimplefp32 = 5.0f0
    argvec = [2.0 3.0 4.0]
    n = [3 6]
    r1 = [125.0 15625.0]
    r2 = [125.0f0 15625.0f0]
    r3 = [8.0 27.0 64.0; 64.0 729.0 4096.0]
    for i in 1:2
        res0 = Swalbe.power_broad(argsimple, n[i])
        @test res0 == r1[i]
        res1 = Swalbe.power_broad(argsimplefp32, n[i])
        @test res1 == r2[i]
        res2 = Swalbe.power_broad.(argvec, n[i])
        for j in 1:3
            @test res2[j] .== r3[i,j]
        end
    end
end