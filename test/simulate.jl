@testset "Simulations" begin
    T = Float64
    @testset "Flat Interface" begin
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=200, tdump=100)
        h = Swalbe.run_flat(sys, "CPU", verbos=false)
        @test all(h .== 1.0)
        @test sum(h) == 25*25
        h = Swalbe.run_flat(sys, "CPU", verbos=true)
        # One dim version
        sys = Swalbe.SysConst_1D(L=25, Tmax=200, tdump=100)
        h = Swalbe.run_flat(sys, verbos=false)
        @test all(h .== 1.0)
        @test sum(h) == 25
        h = Swalbe.run_flat(sys, verbos=true)
    end
    @testset "Random Interface" begin
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=10000, tdump=5000)
        h = Swalbe.run_random(sys, "CPU", ϵ=0.1, verbos=false)
        difference = maximum(h) - minimum(h)
        @test difference < 0.02
        sys = Swalbe.SysConst(Lx=25, Ly=25, Tmax=100, tdump=500)
        h = Swalbe.run_random(sys, "CPU", ϵ=0.1, verbos=true)
        # One dim version
        sys = Swalbe.SysConst_1D(L=25, Tmax=10000, tdump=5000)
        h = Swalbe.run_random(sys, ϵ=0.1, verbos=false)
        difference = maximum(h) - minimum(h)
        @test difference < 0.02
        sys = Swalbe.SysConst_1D(L=25, Tmax=100, tdump=500)
        h = Swalbe.run_random(sys, ϵ=0.1, verbos=true)
    end
    @testset "Rayleigh Taylor" begin
        sys = Swalbe.SysConst(Lx=100, Ly=100, Tmax=1000, tdump=500, g=-0.002)
        _, diff = Swalbe.run_rayleightaylor(sys, "CPU", kx=4, ky=5, ϵ=0.01, verbos=true)
        @test diff[1] < diff[end]
        # TODO: Write a test to check the growthrate 
        # One dim case
        sys = Swalbe.SysConst_1D(L = 512, Tmax=1000, tdump=500, g=-0.002)
        _, diff = Swalbe.run_rayleightaylor(sys, k=4, ϵ=0.01, verbos=true)
        @test diff[1] < diff[end]
    end
    @testset "Relaxing droplet" begin
        rads = 35
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=10000, δ=3.0)
        h, A = Swalbe.run_dropletrelax(sys, "CPU", radius=rads, verbos=false)
        # Initial droplet volume
        vol = π/3 * rads^3 * (2+cospi(1/6)) * (1-cospi(1/6))^2
        R1 = cbrt((rads^3*(2+cospi(1/6))*(1-cospi(1/6))^2)/((2+cospi(1/9))*(1-cospi(1/9))^2))
        r1 = sinpi(1/9)*R1
        # Identify the droplet, below this height is precursor layer
        drop1d = findall(h[75,:] .> 0.055)
        # Get the radius and height
        droprad = length(drop1d)/2
        droph = maximum(h)
        # Compute the droplets shape after 10k time steps
        vnum = 1/6*π*droph*(3*droprad^2 + droph^2)
        # Test that the droplet volume has not changed too much
        @test vol ≈ vnum atol = vol/100*10
        @test r1 ≈ droprad atol = r1/100*10
        @test A[1] < A[end]
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=100, δ=3.0, tdump=50)
        h, A = Swalbe.run_dropletrelax(sys, "CPU", radius=rads, verbos=true)
        # One dim case
        thresh = 0.05026
        sys = Swalbe.SysConst_1D(L = 150, Tmax=10000, δ=3.0)
        h, A = Swalbe.run_dropletrelax(sys, radius=rads, verbos=false)
        # Initial droplet volume
        vol = π/3 * rads^3 * (2+cospi(1/6)) * (1-cospi(1/6))^2
        R1 = cbrt((rads^3*(2+cospi(1/6))*(1-cospi(1/6))^2)/((2+cospi(1/9))*(1-cospi(1/9))^2))
        r1 = sinpi(1/9)*R1
        # Identify the droplet, below this height is precursor layer
        drop1d = findall(h .> 0.055)
        # Drop area (summation of height)
        drop_ar = zeros(sys.L)
        for i in 1:sys.L
            if h[i] > thresh
                drop_ar[i] = h[i]
            end
        end
        # Get the radius and height
        droprad = length(drop1d)/2
        droph = maximum(h)
        # Compute the droplets shape after 10k time steps
        vnum = 1/6*π*droph*(3*droprad^2 + droph^2)
        sph_r = droph/2 + length(findall(h .> thresh))^2/(8*droph)
        ar = sph_r^2/2*(π/9 - sinpi(1/9))
        # Test that the droplet volume has not changed too much
        # TODO fix these two tests
        # @test vol ≈ vnum atol = vol/100*10 # This breaks in 1d I have to investigate at some point in time
        # @test ar ≈ sum(drop_ar) atol = ar/100*10
        @test r1 ≈ droprad atol = r1/100*10
        @test A[1] < A[end]
        sys = Swalbe.SysConst_1D(L = 200, Tmax=100, δ=3.0, tdump=50)
        h, A = Swalbe.run_dropletrelax(sys, radius=35, verbos=true)
    end
    @testset "Patterned droplet" begin
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=10000, δ=3.0)
        h = Swalbe.run_dropletpatterned(sys, "CPU", radius=35, θₛ = fill(1/9,sys.Lx, sys.Ly), verbos=false)
        # Initial droplet volume
        vol = π/3 * 35^3 * (2+cospi(1/6)) * (1-cospi(1/6))^2
        R1 = cbrt((35^3*(2+cospi(1/6))*(1-cospi(1/6))^2)/((2+cospi(1/9))*(1-cospi(1/9))^2))
        r1 = sinpi(1/9)*R1
        # Identify the droplet, below this height is precursor layer
        drop1d = findall(h[75,:] .> 0.055)
        # Get the radius and height
        droprad = length(drop1d)/2
        droph = maximum(h)
        # Compute the droplets shape after 10k time steps
        vnum = 1/6*π*droph*(3*droprad^2 + droph^2)
        # Test that the droplet volume has not changed too much
        @test vol ≈ vnum atol = vol/100*10
        @test r1 ≈ droprad atol = r1/100*10
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=100, tdump=50, δ=3.0)
        h = Swalbe.run_dropletpatterned(sys, "CPU", radius=35, θₛ = fill(1/9,sys.Lx, sys.Ly))
        # One dim case
        thresh = 0.05026
        sys = Swalbe.SysConst_1D(L=150, Tmax=10000, δ=3.0)
        h = Swalbe.run_dropletpatterned(sys, radius=35, θₛ = fill(1/9,sys.L), verbos=false)
        # Initial droplet volume
        vol = π/3 * 35^3 * (2+cospi(1/6)) * (1-cospi(1/6))^2
        R1 = cbrt((35^3*(2+cospi(1/6))*(1-cospi(1/6))^2)/((2+cospi(1/9))*(1-cospi(1/9))^2))
        r1 = sinpi(1/9)*R1
        # Identify the droplet, below this height is precursor layer
        drop1d = findall(h .> 0.055)
        # Drop area (summation of height)
        drop_ar = zeros(sys.L)
        for i in 1:sys.L
            if h[i] > thresh
                drop_ar[i] = h[i]
            end
        end
        # Get the radius and height
        droprad = length(drop1d)/2
        droph = maximum(h)
        # Compute the droplets shape after 10k time steps
        vnum = 1/6*π*droph*(3*droprad^2 + droph^2)
        sph_r = droph/2 + length(findall(h .> thresh))^2/(8*droph)
        ar = sph_r^2/2*(π/9 - sinpi(1/9))
        # TODO fix me
        # @test ar ≈ sum(drop_ar) atol = ar/100*10
        @test r1 ≈ droprad atol = r1/100*10
        sys = Swalbe.SysConst_1D(L=150, Tmax=100, tdump=50, δ=3.0)
        h = Swalbe.run_dropletpatterned(sys, radius=35, θₛ = fill(1/9,sys.L))
    end

    @testset "Sliding droplet" begin
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=5000, δ=2.0)
        # This will return the hydrodynamic moments: h, u, v
        mom = Swalbe.run_dropletforced(sys, "CPU", radius=35, fx=1e-4, verbos=false)
        # Now check that the center has moved 
        @test findmax(mom[1])[2][1] ≠ 75
        @test findmax(mom[1])[2][2] == 75
        # Further test that the velocity is not too large
        @test all(mom[2] .< 0.1) 
        @test all(mom[3] .< 0.1) 
        sys = Swalbe.SysConst(Lx=150, Ly=150, Tmax=100, tdump=50, δ=3.0)
        h = Swalbe.run_dropletforced(sys, "CPU", radius=35, fx=1e-4)
        # One dim case
        sys = Swalbe.SysConst_1D(L=150, Tmax=5000, δ=2.0)
        # This will return the hydrodynamic moments: h, u, v
        mom = Swalbe.run_dropletforced(sys, radius=35, f=1e-4, verbos=false)
        # Now check that the center has moved 
        @test findmax(mom[1])[2] ≠ Int(sys.L/2)
        # Further test that the velocity is not too large
        @test all(mom[2] .< 0.1) 
        sys = Swalbe.SysConst_1D(L=150, Tmax=100, tdump=50, δ=3.0)
        h = Swalbe.run_dropletforced(sys, radius=35, f=1e-4)
    end
end