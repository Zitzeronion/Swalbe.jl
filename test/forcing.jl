@testset "Forcings" begin
    fx = zeros(5,5)
    fy = zeros(5,5)
    f1 = zeros(30)
    # For structs
    sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs())
    state = Swalbe.Sys(sys, "CPU")
    state2 = Swalbe.Sys(sys, "CPU", kind="thermal")
    # One dim model   
    sys1D = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs())
    st1 = Swalbe.Sys(sys1D)
    st2 = Swalbe.Sys(sys1D, kind="thermal")
    @testset "Slippage" begin
        # No velocities
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), zeros(5,5), 1.0, 1/6)
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state2, sys)
        for i in [(fx, fy), (state.slipx, state.slipy), (state2.basestate.slipx, state2.basestate.slipy)]
            @test all(i[1] .== 0.0)
            @test all(i[2] .== 0.0)
        end
        
        # Velocity in x
        state.velx .= 0.1
        state2.basestate.velx .= 0.1
        Swalbe.slippage!(fx, fy, ones(5,5), fill(0.1,5,5), zeros(5,5), 1.0, 1/6)
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state2, sys)
        for i in [(fx, fy), (state.slipx, state.slipy), (state2.basestate.slipx, state2.basestate.slipy)]
            @test all(isapprox.(i[1], 0.1/11; atol=1e-10))
            @test all(i[2] .== 0.0)
        end
        
        # Velocity in y
        for i in [state, state2.basestate]
            i.velx .= 0.0
            i.vely .= 0.1
        end
        Swalbe.slippage!(fx, fy, ones(5,5), zeros(5,5), fill(0.1,5,5), 1.0, 1/6)
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state2, sys)
        for i in [(fx, fy), (state.slipx, state.slipy), (state2.basestate.slipx, state2.basestate.slipy)]
            @test all(isapprox.(i[2], 0.1/11; atol=1e-10))
            @test all(i[1] .== 0.0)
        end
        
        # Velocity
        for i in [state, state2.basestate]
            i.velx .= -0.1
            i.vely .= 0.1
        end
        Swalbe.slippage!(fx, fy, ones(5,5), fill(-0.1,5,5), fill(0.1,5,5), 1.0, 1/6)
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state2, sys)
        for i in [(fx, fy), (state.slipx, state.slipy), (state2.basestate.slipx, state2.basestate.slipy)]
            @test all(isapprox.(i[1], -0.1/11; atol=1e-10))
            @test all(isapprox.(i[2], 0.1/11; atol=1e-10))
        end
        
        # No slip
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(δ=0.0))
        Swalbe.slippage!(fx, fy, ones(5,5), fill(-0.1,5,5), fill(0.1,5,5), 0.0, 1/6)
        Swalbe.slippage!(state, sys)
        Swalbe.slippage!(state2, sys)
        for i in [(fx, fy), (state.slipx, state.slipy), (state2.basestate.slipx, state2.basestate.slipy)]
            @test all(isapprox.(i[1], -0.1/2; atol=1e-10))
            @test all(isapprox.(i[2], 0.1/2; atol=1e-10))
        end
    end

    @testset "Slippage 1D" begin
        # No velocities
        Swalbe.slippage!(f1, ones(30), zeros(30), 1.0, 1/6)
        Swalbe.slippage!(st1, sys1D)
        Swalbe.slippage!(st2, sys1D)
        for i in [f1, st1.slip, st2.basestate.slip]
            @test all(i .== 0.0)
        end
        # Velocity in x
        for i in [st1, st2.basestate]
            i.vel .= 0.1
        end
        Swalbe.slippage!(f1, ones(30), fill(0.1,30), 1.0, 1/6)
        Swalbe.slippage!(st1, sys1D)
        Swalbe.slippage!(st2, sys1D)
        for i in [f1, st1.slip, st2.basestate.slip]
            @test all(isapprox.(i, 0.1/11; atol=1e-10))
        end            
        
        # No slip
        sys1D = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(δ=0.0))
        for i in [st1, st2.basestate]
            i.vel .= -0.1
        end
        Swalbe.slippage!(f1, ones(30), fill(-0.1,30), 0.0, 1/6)
        Swalbe.slippage!(st1, sys1D)
        Swalbe.slippage!(st2, sys1D)
        for i in [f1, st1.slip, st2.basestate.slip]
            @test all(isapprox.(i, -0.1/2; atol=1e-10))
        end
        
    end
    
    @testset "Pressure gradient" begin
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs())
        sys1 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs())
        state = Swalbe.Sys(sys, "CPU")    
        state2 = Swalbe.Sys(sys, "CPU", kind="thermal")  
        st = Swalbe.Sys(sys1)    
        st2 = Swalbe.Sys(sys1, kind="thermal")     
        state.pressure .= reshape(collect(1.0:25),5,5)
        state2.basestate.pressure .= reshape(collect(1.0:25),5,5)
        st.pressure .= collect(1.0:30)
        st2.basestate.pressure .= collect(1.0:30)
        Swalbe.h∇p!(state)
        Swalbe.h∇p!(state2)
        Swalbe.h∇p!(st)
        Swalbe.h∇p!(st2)
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

        sol = ones(30)
        sol[1] = sol[end] = -14
        for i in [(state.h∇px, state.h∇py), (state2.basestate.h∇px, state2.basestate.h∇py)]
            @test all(i[1] .== solx)
            @test all(i[2] .== soly)
        end
        for i in [st.h∇p, st2.basestate.h∇p]
            @test all(i .== sol)
        end
    end
    
    @testset "Thermal" begin
        f1 = ones(50,50)
        f2 = ones(50,50)
        f1D = ones(100000)
        
        sys = Swalbe.SysConst(Lx=50, Ly=50, param=Swalbe.Taumucs(kbt=0.01))
        sys1 = Swalbe.SysConst_1D(L=100000, param=Swalbe.Taumucs(kbt=0.01))
        state = Swalbe.Sys(sys, "CPU", kind="thermal")    
        st1 = Swalbe.Sys(sys1, kind="thermal")    
        
        for kb in [0.01, 0.1]
            sys = Swalbe.SysConst(Lx=50, Ly=50, param=Swalbe.Taumucs(kbt=kb))
            sys1 = Swalbe.SysConst_1D(L=100000, param=Swalbe.Taumucs(kbt=kb))
            vartest = 2*kb/11
            Swalbe.thermal!(f1, f2, ones(50,50), kb, 1/6, 1.0)
            Swalbe.thermal!(f1D, ones(100000), kb, 1/6, 1.0)
            Swalbe.thermal!(state, sys)
            Swalbe.thermal!(st1, sys1)
            for i in enumerate([f1, f2, state.kbtx, state.kbty, st1.kbt])
                @test mean(i[2]) ≈ 0.0 atol=1e-2
                @test var(i[2]) ≈ vartest atol=vartest/10
            end
        end
    end

    @testset "inclination" begin
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs())
        sys1 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs())
        state = Swalbe.Sys(sys, "CPU")    
        st1 = Swalbe.Sys(sys1) 
        state2 = Swalbe.Sys(sys, "CPU", kind="thermal")    
        st2 = Swalbe.Sys(sys1, kind="thermal")
        sols = Dict(0=>0.05, 1=>0.1*(0.5 + 0.5 * tanh(1.0)))
        for t in [0, 1]
            Swalbe.inclination!([0.1, 0.1], state, t=t, tstart=0, tsmooth=1)
            Swalbe.inclination!([0.1, 0.1], state2, t=t, tstart=0, tsmooth=1)
            Swalbe.inclination!(0.1, st1, t=t, tstart=0, tsmooth=1)
            Swalbe.inclination!(0.1, st2, t=t, tstart=0, tsmooth=1)
            for i in [state.Fx, state.Fy, state2.basestate.Fx, state2.basestate.Fy, st1.F, st2.basestate.F]
                @test all(i .== sols[t])
                i .= 0.0
            end
        end
    end
  
    @testset "Surface tension gradient" begin
        sys1 = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs())
        state = Swalbe.Sys(sys1, kind="gamma") 
        state.γ .= collect(1.0:30) 
        Swalbe.∇γ!(state)
        sol = fill(3/2, 30)
        sol[1] = sol[end] = -21.0
        @test all(state.∇γ .== sol)
    end

    @testset "rho update" begin
        @testset "Constant fields" begin
            rho = ones(25)
            height = ones(25)
            output = zeros(25)
            Swalbe.update_rho!(rho, output, height, zeros(25,2), zeros(25,4))
            @test all(rho .== 1)
        end
    end
  
    @testset "view four" begin
        dummy = reshape(collect(1:20),5,4)
        d1, d2, d3, d4 = Swalbe.view_four(dummy)
        @test all(d1 .== dummy[:,1])
        @test all(d2 .== dummy[:,2])
        @test all(d3 .== dummy[:,3])
        @test all(d4 .== dummy[:,4])
    end
end