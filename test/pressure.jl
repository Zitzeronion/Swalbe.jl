@testset "Capillary pressure" begin
    # Struct
    sys = Swalbe.SysConst(Lx=5, Ly=5, n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.1)
    state = Swalbe.Sys(sys, "CPU")
    state.height .= reshape(collect(1.0:25),5,5)
    # With arguments struct
    f = reshape(collect(1.0:25),5,5)
    f_float = reshape(collect(1.0f0:25.0f0),5,5)
    res = zeros(5,5)
    @testset "No contact angle" begin
        Swalbe.filmpressure!(res, f, 1.0, 0.0, 3, 2, 0.1, 0.1)
        Swalbe.filmpressure!(state, sys, 0.0)
        sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0;
               ]
        for i in eachindex(sol)
            @test res[i] .≈ sol[i] atol=1e-10
            @test state.pressure[i] .≈ sol[i] atol=1e-10
        end
        Swalbe.filmpressure!(res, f, 0.0)
        for i in eachindex(sol)
            @test res[i] .≈ 0.01 * sol[i] atol=1e-10
        end
    end
    @testset "No height gradient" begin
        nograd = ones(5,5)
        state.height .= 1.0
        sys2 = Swalbe.SysConst(Lx=5, Ly=5, n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0)
        Swalbe.filmpressure!(res, nograd, 1.0, 1/2, 3, 2, 0.1, 0.0)
        Swalbe.filmpressure!(state, sys2, 1/2)
        for i in eachindex(res)
            @test res[i] .≈ -2(0.1^2-0.1) atol=1e-10
            @test state.pressure[i] .≈ -2(0.1^2-0.1) atol=1e-10
        end
    end
    @testset "Gradient and contact angle" begin
        state.height .= reshape(collect(1.0:25),5,5)
        sys2 = Swalbe.SysConst(Lx=5, Ly=5, n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0)
        Swalbe.filmpressure!(state, sys2, 1/2)
        Swalbe.filmpressure!(res, f, 1.0, 1/2, 3, 2, 0.1, 0.0)
        sol = [30.0 5.0 5.0 5.0 -20;
               25.0 0.0 0.0 0.0 -25.0;
               25.0 0.0 0.0 0.0 -25.0;
               25.0 0.0 0.0 0.0 -25.0;
               20.0 -5.0 -5.0 -5.0 -30.0]
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            @test res[i] .≈ -1(sol[i] + 20((0.1/f[i])^3-(0.1/f[i])^2)) atol=1e-10
            @test state.pressure[i] .≈ -1(sol[i] + 20((0.1/f[i])^3-(0.1/f[i])^2)) atol=1e-10
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
            @test res[i] .≈ -1(sol[i] + 20((0.1f0/f_float[i])^3-(0.1f0/f_float[i])^2)) atol=1e-6
        end
    end
    dgrad = zeros(5,5,8)
    @testset "No contact angle circshift!" begin
        state.height .= reshape(collect(1.0:25),5,5)
        Swalbe.filmpressure!(state, sys, 0.0)
        Swalbe.filmpressure!(res, f, dgrad, 1.0, 0.0, 3, 2, 0.1, 0.1)
        sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0;
               ]
        for i in eachindex(sol)
            @test res[i] .≈ sol[i] atol=1e-10
            @test state.pressure[i] .≈ sol[i] atol=1e-10
        end
    end
end

@testset "Capillary pressure 1D" begin
    # Struct
    sys = Swalbe.SysConst_1D(L=30, n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.1)
    state = Swalbe.Sys(sys)
    state.height .= collect(1.0:30)
    sys2 = Swalbe.SysConst_1D(L=30, n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0)
    state2 = Swalbe.Sys(sys2)
    # Without the struct
    f = collect(1.0:30)
    f_float = collect(1.0f0:30.0f0)
    res = zeros(30)
    dummy = zeros(30,3)
    @testset "No contact angle" begin
        Swalbe.filmpressure!(res, f, dummy, 1.0, 0.0, 3, 2, 0.1, 0.1)
        Swalbe.filmpressure!(state, sys, 0.0)
        # println("My result: $res")
        sol = zeros(30)
        sol[1] = 30
        sol[end] = -30
        @test all(res .== -sol)
        @test all(state.pressure .== -sol)
    end
    @testset "No height gradient" begin
        nograd = ones(30)
        state2.height .= 1.0
        Swalbe.filmpressure!(res, nograd, dummy, 1.0, 1/2, 3, 2, 0.1, 0.0)
        Swalbe.filmpressure!(state2, sys2, 1/2)
        @test res[15] ≈ -2(0.1^2-0.1) atol=1e-10
        @test state2.pressure[15] ≈ -2(0.1^2-0.1) atol=1e-10
        
    end
    @testset "Gradient and contact angle" begin
        Swalbe.filmpressure!(res, f, dummy, 1.0, 1/2, 3, 2, 0.1, 0.0)
        Swalbe.filmpressure!(state, sys2, 1/2)
        sol = zeros(30)
        sol[1] = 30
        sol[end] = -30
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            @test res[i] .≈ -1(sol[i] + 20((0.1/f[i])^3-(0.1/f[i])^2)) atol=1e-10
            @test state.pressure[i] .≈ -1(sol[i] + 20((0.1/f[i])^3-(0.1/f[i])^2)) atol=1e-10
        end
    end
end

@testset "Capillary pressure 1D with Gamma" begin
    f = collect(1.0:30)
    res = zeros(30)
    rho = fill(0.1, 30)
    dummy = zeros(30,2)
    @testset "No Gamma no theta" begin
        Swalbe.filmpressure!(res, f, dummy, rho, 1.0, 0.0, 3, 2, 0.1, 0.1)
        # println("My result: $res")
        sol = zeros(30)
        sol[1] = 30
        sol[end] = -30
        @test all(res .== -sol)
    end
    @testset "No Gamma no laplacian" begin
        nograd = ones(30)
        Swalbe.filmpressure!(res, nograd, dummy, rho, 1.0, 1/2, 3, 2, 0.1, 0.0)
        for i in eachindex(res)
            @test res[i] .≈ -2(0.1^2-0.1) atol=1e-10
        end
    end
    @testset "Gamma but rho=0 no laplacian" begin
        nograd = ones(30)
        Swalbe.filmpressure!(res, nograd, dummy, zeros(30), 1.0, 1/2, 3, 2, 0.1, 0.0)
        for i in eachindex(res)
            @test res[i] .≈ -2(0.1^2-0.1) atol=1e-10
        end
    end
    @testset "Gamma, rho, gradient and contact angle" begin
        rho = fill(0.1, 30)
        Swalbe.filmpressure!(res, f, dummy, rho, 1.0, 1/2, 3, 2, 0.1, 0.0, Gamma=0.1)
        sol = zeros(30)
        sol[1] = 30
        sol[end] = -30
        for i in eachindex(res)
            # Now the result comprisses of two components the disjoining potential and the laplace term.
            # The prefactor (1 + 0.01) -> (gamma + Gamma * rho)
            @test res[i] .≈ -(1 + 0.01) * (sol[i] + 20((0.1/f[i])^3-(0.1/f[i])^2)) atol=1e-10
        end
    end
end

@testset "Power function" begin
    arg = 2
    argsimple = 5.0
    argsimplefp32 = 5.0f0
    argvec = [2.0 3.0 4.0]
    n = [3 6]
    r = [8 64]
    r1 = [125.0 15625.0]
    r2 = [125.0f0 15625.0f0]
    r3 = [8.0 27.0 64.0; 64.0 729.0 4096.0]
    for i in 1:2
        res0 = Swalbe.power_broad(arg, n[i])
        @test res0 == r[i] 
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