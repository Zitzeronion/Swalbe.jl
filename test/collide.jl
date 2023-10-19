@testset "Collision and Stream" begin
    function set_dist!(a, b, c)
        a .= 1.0
        b .= 1.0
        c .= 1.0
        a[1,1,:] .= 2.0
        return nothing
    end
    # Allocations
    feq = ones(5,5,9)
    ftemp = ones(5,5,9)
    fout = ones(5,5,9)
    feq[1,1,:] .= 2.0
    # States
    sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(τ=0.75))
    state = Swalbe.Sys(sys, "CPU")
    state2 = Swalbe.Sys(sys, "CPU", kind="thermal")
    for i in [state, state2.basestate]
        i.feq .= 1.0
        i.ftemp .= 1.2
        i.fout .= 1.0
        i.feq[1,1,:] .= 2.0
    end

    onebytau = 1.0/0.75
    omega = 1.0 - 1.0/0.75
    @testset "Dummy dists τ=1 no forces" begin
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(τ=1.0))
        Swalbe.BGKandStream!(fout, feq, ftemp, zeros(5,5), zeros(5,5), 1.0)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.BGKandStream!(state2, sys)
        for i in [(fout, feq), (state.fout, state.feq), (state2.basestate.fout, state2.basestate.feq)]
            @test all(i[1][:,:,1] .== i[2][:,:,1])
            @test all(i[1][:,:,2] .== circshift(i[2][:,:,1],(1,0)))
            @test all(i[1][:,:,3] .== circshift(i[2][:,:,1],(0,1)))
            @test all(i[1][:,:,4] .== circshift(i[2][:,:,1],(-1,0)))
            @test all(i[1][:,:,5] .== circshift(i[2][:,:,1],(0,-1)))
            @test all(i[1][:,:,6] .== circshift(i[2][:,:,1],(1,1)))
            @test all(i[1][:,:,7] .== circshift(i[2][:,:,1],(-1,1)))
            @test all(i[1][:,:,8] .== circshift(i[2][:,:,1],(-1,-1)))
            @test all(i[1][:,:,9] .== circshift(i[2][:,:,1],(1,-1)))
        end
    end
    @testset "Dummy dists τ=1 with forces" begin
        for i in [(feq, ftemp, fout), (state.feq, state.ftemp, state.fout), 
            (state2.basestate.feq, state2.basestate.ftemp, state2.basestate.fout)]
            set_dist!(i[1], i[2], i[3])
        end
        feq .= 1.0
        for i in [state, state2.basestate]
            i.feq .= 1.0
            i.ftemp .= 1.2
            i.Fx .= 0.1
            i.Fy .= -0.1
        end
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(τ=1.0))
        Swalbe.BGKandStream!(state, sys)
        Swalbe.BGKandStream!(state2, sys)
        Swalbe.BGKandStream!(fout, feq, ftemp, fill(0.1,5,5), fill(-0.1,5,5), 1.0)
        for i in [(fout, feq), (state.fout, state.feq), (state2.basestate.fout, state2.basestate.feq)]
            @test all(i[1][:,:,1] .== i[2][:,:,1])
            @test all(i[1][:,:,2] .== circshift(i[2][:,:,1] .+ 1/30,(1,0)))
            @test all(i[1][:,:,3] .== circshift(i[2][:,:,1] .- 1/30,(0,1)))
            @test all(i[1][:,:,4] .== circshift(i[2][:,:,1] .- 1/30,(-1,0)))
            @test all(i[1][:,:,5] .== circshift(i[2][:,:,1] .+ 1/30,(0,-1)))
            @test all(i[1][:,:,6] .== circshift(i[2][:,:,1],(1,1)))
            @test all(i[1][:,:,7] .== circshift(i[2][:,:,1] .- 1/24 .* 0.2 ,(-1,1)))
            @test all(i[1][:,:,8] .== circshift(i[2][:,:,1],(-1,-1)))
            @test all(i[1][:,:,9] .== circshift(i[2][:,:,1] .+ 1/24 .* 0.2,(1,-1)))
        end
    end
    @testset "Dummy dists τ=0.75 no forces" begin
        
        for i in [(feq, ftemp, fout), (state.feq, state.ftemp, state.fout), 
                  (state2.basestate.feq, state2.basestate.ftemp, state2.basestate.fout)]
            set_dist!(i[1], i[2], i[3])
        end
        for i in [state, state2.basestate]
            i.Fx .= 0.0
            i.Fy .= 0.0
        end
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(τ=0.75))
        Swalbe.BGKandStream!(fout, feq, ftemp, zeros(5,5), zeros(5,5), 0.75)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.BGKandStream!(state2, sys)
        for i in [(fout, feq), (state.fout, state.feq), (state2.basestate.fout, state2.basestate.feq)]
            @test all(i[1][:,:,1] .== omega .* 1.0 .+ onebytau .* i[2][:,:,1])
            @test all(i[1][:,:,2] .≈ circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(1,0)))
            @test all(i[1][:,:,3] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(0,1)))
            @test all(i[1][:,:,4] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(-1,0)))
            @test all(i[1][:,:,5] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(0,-1)))
            @test all(i[1][:,:,6] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(1,1)))
            @test all(i[1][:,:,7] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(-1,1)))
            @test all(i[1][:,:,8] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(-1,-1)))
            @test all(i[1][:,:,9] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(1,-1)))
        end
    end
    
    @testset "Dummy dists τ=0.75 with forces" begin
        for i in [(feq, ftemp, fout), (state.feq, state.ftemp, state.fout), 
            (state2.basestate.feq, state2.basestate.ftemp, state2.basestate.fout)]
            set_dist!(i[1], i[2], i[3])
        end
        for i in [state, state2.basestate]
            i.Fx .= 0.1
            i.Fy .= -0.1
        end
        sys = Swalbe.SysConst(Lx=5, Ly=5, param=Swalbe.Taumucs(τ=0.75))
        Swalbe.BGKandStream!(fout, feq, ftemp, fill(0.1,5,5), fill(-0.1,5,5), 0.75)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.BGKandStream!(state2, sys)
        for i in [(fout, feq), (state.fout, state.feq), (state2.basestate.fout, state2.basestate.feq)]
            @test all(i[1][:,:,1] .== omega .* 1.0 .+ onebytau .* i[2][:,:,1])
            @test all(i[1][:,:,2] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1].+ 1/30,(1,0)))
            @test all(i[1][:,:,3] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1].- 1/30,(0,1)))
            @test all(i[1][:,:,4] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1].- 1/30,(-1,0)))
            @test all(i[1][:,:,5] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1].+ 1/30,(0,-1)))
            @test all(i[1][:,:,6] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(1,1)))
            @test all(i[1][:,:,7] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1] .- 1/24 .* 0.2,(-1,1)))
            @test all(i[1][:,:,8] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1],(-1,-1)))
            @test all(i[1][:,:,9] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,:,1] .+ 1/24 .* 0.2,(1,-1)))
        end
    end

    feq = ones(30,3)
    ftemp = ones(30,3)
    fout = ones(30,3)
    feq[1,:] .= 2.0
    sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(τ=0.75))
    st1d = Swalbe.Sys(sys)
    st1d2 = Swalbe.Sys(sys, kind="thermal")
    # Using a wall
    obst = zeros(30 )
    obst[1]=1
    obst[end]=1
    sysBC = Swalbe.SysConstWithBound_1D{Float64}(obs=obst, L=30, param=Swalbe.Taumucs(τ=0.75))
    Swalbe.obslist!(sysBC)
    st1dBC = Swalbe.Sys(sysBC, kind="gamma_bound")
    onebytau = 1.0/0.75
    omega = 1.0 - 1.0/0.75
    @testset "Dummy dists τ=1 no forces 1D" begin
        sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(τ=1.0))
        sysBC = Swalbe.SysConstWithBound_1D{Float64}(obs=obst, L=30, param=Swalbe.Taumucs(τ=1.0))
        Swalbe.BGKandStream!(fout, feq, ftemp, zeros(30), 1.0)
        Swalbe.BGKandStream!(st1d, sys)
        Swalbe.BGKandStream!(st1d2, sys)
        Swalbe.BGKandStream!(st1dBC, sysBC)
        for i in [(fout, feq), (st1d.fout, st1d.feq), (st1d2.basestate.fout, st1d2.basestate.feq), 
            (st1dBC.basestate.fout, st1dBC.basestate.feq)]
            @test all(i[1][:,1] .== i[2][:,1])
            @test all(i[1][:,2] .== circshift(i[2][:,1],1))
            @test all(i[1][:,3] .== circshift(i[2][:,1],-1))
        end
    end
    @testset "Dummy dists τ=0.75 no forces 1D" begin
        for i in [(feq, ftemp, fout), (st1d.feq, st1d.ftemp, st1d.fout), 
            (st1d2.basestate.feq, st1d2.basestate.ftemp, st1d2.basestate.fout)]
            i[1] .= 1.0
            i[2] .= 1.0
            i[3] .= 1.0
            i[1][1,:] .= 2.0
        end
        sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(τ=0.75))
        Swalbe.BGKandStream!(fout, feq, ftemp, zeros(30), 0.75)
        Swalbe.BGKandStream!(st1d, sys)
        Swalbe.BGKandStream!(st1d2, sys)
        for i in [(fout, feq), (st1d.fout, st1d.feq), (st1d2.basestate.fout, st1d2.basestate.feq)]
            @test all(i[1][:,1] .== omega .* 1.0 .+ onebytau .* i[2][:,1])
            @test all(i[1][:,2] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,1],1))
            @test all(i[1][:,3] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,1],-1))
        end
    end
    
    @testset "Dummy dists τ=1 with forces 1D" begin
        for i in [(feq, ftemp, fout), (st1d.feq, st1d.ftemp, st1d.fout), 
            (st1d2.basestate.feq, st1d2.basestate.ftemp, st1d2.basestate.fout)]
            i[1] .= 1.0
            i[2] .= 1.0
            i[3] .= 1.0
            i[1][1,:] .= 2.0
        end
        for i in [st1d, st1d2.basestate]
            i.F .= 0.1
        end
        sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(τ=1.0))
        Swalbe.BGKandStream!(fout, feq, ftemp, fill(0.1,30), 1.0)
        Swalbe.BGKandStream!(st1d, sys)
        Swalbe.BGKandStream!(st1d2, sys)
        for i in [(fout, feq), (st1d.fout, st1d.feq), (st1d2.basestate.fout, st1d2.basestate.feq)]
            @test all(i[1][:,1] .== i[2][:,1])
            @test all(i[1][:,2] .== circshift(i[2][:,1] .+ 1/20,1))
            @test all(i[1][:,3] .== circshift(i[2][:,1] .- 1/20,-1))
        end
    end
    
    @testset "Dummy dists τ=0.75 with forces 1D" begin
        for i in [(feq, ftemp, fout), (st1d.feq, st1d.ftemp, st1d.fout), 
            (st1d2.basestate.feq, st1d2.basestate.ftemp, st1d2.basestate.fout)]
            i[1] .= 1.0
            i[2] .= 1.0
            i[3] .= 1.0
            i[1][1,:] .= 2.0
        end
        for i in [st1d, st1d2.basestate]
            i.F .= 0.1
        end
        sys = Swalbe.SysConst_1D(L=30, param=Swalbe.Taumucs(τ=0.75))
        Swalbe.BGKandStream!(fout, feq, ftemp, fill(0.1,30), 0.75)
        Swalbe.BGKandStream!(st1d, sys)
        Swalbe.BGKandStream!(st1d2, sys)
        for i in [(fout, feq), (st1d.fout, st1d.feq), (st1d2.basestate.fout, st1d2.basestate.feq)]
            @test all(i[1][:,1] .== omega .* 1.0 .+ onebytau .* i[2][:,1])
            @test all(i[1][:,2] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,1].+ 1/20,1))
            @test all(i[1][:,3] .== circshift(omega .* 1.0 .+ onebytau * i[2][:,1].- 1/20,-1))
        end
    end
end

@testset "viewdists" begin
    f = reshape(collect(1.0:225.0),5,5,9)
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(f)
    allviews = [f0, f1, f2, f3, f4, f5, f6, f7, f8]
    for (index, value) in enumerate(allviews)
        @test all(value .== f[:,:,index])
    end
end

@testset "viewdists_1D" begin
    f = reshape(collect(1.0:30.0),10,3)
    f0, f1, f2 = Swalbe.viewdists_1D(f)
    allviews = [f0, f1, f2]
    for (index, value) in enumerate(allviews)
        @test all(value .== f[:,index])
    end
end