@testset "Simulations" begin
    @testset "Flat Interface" begin
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=200)
        h = Swalbe.run_flat(sys, "CPU", verbos=false)
        @test all(h .== 1.0)
        @test sum(h) == 25*25
    end
    @testset "Random Interface" begin
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=10000)
        h = Swalbe.run_random(sys, "CPU", ϵ=0.1, verbos=false)
        difference = maximum(h) - minimum(h)
        @test difference < 0.02
    end
    @testset "Relaxing droplet" begin
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=10000)
        h = Swalbe.run_dropletrelax(sys, "CPU", radius=35, verbos=false)
        vol = π/3 * 35^3 * (2+cospi(1/6)) * (1-cospi(1/6))^2
        drop1d = findall(h[75,:] .> 0.055)
        dropdiameter = (drop1d[end] - drop1d[1])/2
        droph = maximum(h)
        vnum = 1/6*π*droph*(3*dropdiameter^2 + droph^2)
        @test vol ≈ vnum atol = vol/100*10
    end
end