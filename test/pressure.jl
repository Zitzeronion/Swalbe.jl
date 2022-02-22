@testset "Capillary pressure" begin
    # Struct
    sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.1))
    state = Swalbe.Sys(sys, "CPU")
    state2 = Swalbe.Sys(sys, "CPU", kind="thermal")
    state.height .= reshape(collect(1.0:25),5,5)
    state2.basestate.height .= reshape(collect(1.0:25),5,5)
    # With arguments struct
    f = reshape(collect(1.0:25),5,5)
    res = zeros(5,5)
    dgrad = zeros(5,5,8)
    sol = [-30.0 -5.0 -5.0 -5.0 20;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -25.0 0.0 0.0 0.0 25.0;
               -20.0 5.0 5.0 5.0 30.0;
               ]
    @testset "No contact angle" begin
        Swalbe.filmpressure!(res, f, dgrad, 1.0, 0.0, 3, 2, 0.1, 0.1)
        Swalbe.filmpressure!(state, sys, θ=0.0)
        Swalbe.filmpressure!(state2, sys, θ=0.0)
        

        @test all(isapprox.(res, sol; atol=1e-10))
        @test all(isapprox.(state.pressure, sol; atol=1e-10))
    end
    @testset "No contact angle circshift!" begin
        state.height .= reshape(collect(1.0:25),5,5)
        Swalbe.filmpressure!(state, sys, θ=0.0)
        Swalbe.filmpressure!(res, f, dgrad, 1.0, 0.0, 3, 2, 0.1, 0.1)
        
        @test all(isapprox.(res, sol; atol=1e-10))
        @test all(isapprox.(state.pressure, sol; atol=1e-10))
    end
    @testset "Gradient and contact angle" begin
        sys2 = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0))
        Swalbe.filmpressure!(state, sys2, θ=1/2)
        Swalbe.filmpressure!(state2, sys2, θ=1/2)
        Swalbe.filmpressure!(res, f, dgrad, 1.0, 1/2, 3, 2, 0.1, 0.0)
        for i in [res, state.pressure, state2.basestate.pressure]
            @test all(isapprox.(i, -1 .* (-sol .+ 20 .* ((0.1 ./ f).^3 .- (0.1 ./ f).^2)); atol=1e-10))
        end
    end
    @testset "No height gradient" begin
        nograd = ones(5,5)
        state.height .= 1.0
        sys2 = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0, θ=1/2))
        Swalbe.filmpressure!(res, nograd, dgrad, 1.0, 1/2, 3, 2, 0.1, 0.0)
        Swalbe.filmpressure!(state, sys2, θ=1/2)
        for i in [res, state.pressure]
            @test all(isapprox.(i, -2(0.1^2-0.1); atol=1e-10))
        end
    end
end

@testset "Capillary pressure 1D" begin
    # Struct
    sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.1))
    state = Swalbe.Sys(sys)
    state.height .= collect(1.0:30)
    sys2 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0, θ=1/2))
    state2 = Swalbe.Sys(sys2)
    sys3 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0, θ=0.0))
    state3 = Swalbe.Sys(sys3, kind="gamma")
    state3.γ .= sys3.param.γ
    state3.basestate.height .= collect(1.0:30)
    # Without the struct
    f = collect(1.0:30)
    sol = zeros(30)
    sol[1] = 30
    sol[end] = -30
    res = zeros(30)
    dummy = zeros(30,3)
    @testset "No contact angle" begin
        Swalbe.filmpressure!(res, f, dummy, 1.0, 0.0, 3, 2, 0.1, 0.1)
        Swalbe.filmpressure!(state, sys, θ=0.0)
        Swalbe.filmpressure!(state3, sys3)
        for i in [res, state.pressure, state3.basestate.pressure]
            @test all(i .== -sol)
        end
    end
    @testset "Gradient and contact angle" begin
        Swalbe.filmpressure!(res, f, dummy, 1.0, 1/2, 3, 2, 0.1, 0.0)
        Swalbe.filmpressure!(state, sys2, θ=1/2)
        sys3 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0, θ=1/2))
        Swalbe.filmpressure!(state3, sys3)
        for i in [res, state.pressure, state3.basestate.pressure]
            @test all(isapprox.(i, -1 .* (sol .+ 20 .* ((0.1 ./ f).^3 .- (0.1 ./ f).^2)); atol=1e-10))
        end
    end
    @testset "No height gradient" begin
        nograd = ones(30)
        state2.height .= 1.0
        state3.basestate.height .= 1.0
        Swalbe.filmpressure!(res, nograd, dummy, 1.0, 1/2, 3, 2, 0.1, 0.0)
        Swalbe.filmpressure!(state2, sys2, θ=1/2)
        sys3 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(n=3, m=2, γ=1.0, hmin=0.1, hcrit=0.0, θ=1/2))
        Swalbe.filmpressure!(state3, sys3)
        for i in [res, state2.pressure, state3.basestate.pressure]
            @test all(isapprox.(i, -2(0.1^2-0.1); atol=1e-10))
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
        @test all(isapprox.(res, -2(0.1^2-0.1); atol=1e-10))
    end
    @testset "Gamma but rho=0 no laplacian" begin
        nograd = ones(30)
        Swalbe.filmpressure!(res, nograd, dummy, zeros(30), 1.0, 1/2, 3, 2, 0.1, 0.0)
        @test all(isapprox.(res, -2(0.1^2-0.1); atol=1e-10))
    end
    @testset "Gamma, rho, gradient and contact angle" begin
        rho = fill(0.1, 30)
        Swalbe.filmpressure!(res, f, dummy, rho, 1.0, 1/2, 3, 2, 0.1, 0.0, Gamma=0.1)
        sol = zeros(30)
        sol[1] = 30
        sol[end] = -30
        @test all(isapprox.(res, -(1 + 0.01) .* (sol .+ 20 .* ((0.1 ./ f).^3 .- (0.1 ./ f).^2)); atol=1e-10))
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