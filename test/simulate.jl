@testset "Simulations" begin
    @testset "Flat Interface" begin
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=200)
        h = Swalbe.run_flat(sys, "CPU")
        @test all(h .== 1.0)
        @test sum(h) == 25*25
    end
    @testset "Random Interface" begin
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=10000)
        h = Swalbe.run_random(sys, "CPU", ϵ=0.1, verbos=false)
        difference = maximum(h) - minimum(h)
        @test difference < 0.02
    end
end