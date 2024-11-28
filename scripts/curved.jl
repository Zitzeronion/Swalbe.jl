using Pkg; Pkg.activate(".")
using Swalbe, DelimitedFiles, Dates, LaTeXStrings
using Plots; Plots.gr(); ENV["GKSwstype"] = "100"
function simulation(;
	Tmax=Int(1e5),                       # Simulation time
        tdump=Int(floor(Tmax/100)),     # dump time
        h= 1,                           # initial film height
        eps=0.001,                      # initial noise
        gamma=1/6^2,                    # surface tension
        data = "data",                  # folder for saving data
        L=2^8,                         # system size
        theta=1/9,                      # contact angle
        n=9, m=3, hmin=0.2,             # parameters for disjoining pressure
        delta=1.0,                      # slip length, in combination with heigher hmin i
        omega=3,                        # number of periods of the underlying substrate
	plot=true			# do plots
)
        # set up system
	sys=Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=Tmax, tdump=tdump, γ = gamma, δ=delta,θ =theta, n=n,m=m,hmin=hmin))
        try
                mkdir(data)
                mkdir("$(data)/figs")
        catch
        end
        state = Swalbe.Sys(sys, kind="curved")
        # Initial data
        state.height.= h .+ eps .* rand(sys.L)
        # set up system profile
        state.substrate .=[1/4*sin(i*omega*2*pi/sys.L) for i in 1:sys.L]
        # We need the gradient of the substrate profile
        Swalbe.∇f!(state.grad_substrate, state.substrate, state.dgrad)
        println("Starting Lattice Boltzmann time loop for active thin film")
        i=0
        for t in 0:sys.param.Tmax
                # save data
                if t % sys.param.tdump == 0
			println("Time step t=$t, mass=$(sum(state.height)), diff=$(maximum(state.height)-minimum(state.height))")
                        writedlm("$(data)/$(Swalbe.to_4_digits(i))_height.csv", state.height)
			if plot
				Plots.plot(state.height .+ state.substrate, label=L"h", c=:cyan, fillrange=minimum(state.substrate))
				Plots.plot!(state.substrate, label=L"s", c=:brown, fillrange=minimum(state.substrate))
				Plots.savefig("$(data)/figs/$(Swalbe.to_4_digits(i)).png")
			end
                        i+=1
                end
                # system update, using the pressure for curved substrates
                Swalbe.filmpressure_curved!(state, sys)
                Swalbe.h∇p!(state)
                Swalbe.slippage!(state, sys)
                state.F .= -state.h∇p .- state.slip
                Swalbe.equilibrium!(state)
                Swalbe.BGKandStream!(state, sys)
                Swalbe.moments!(state)
        end
    println("Finished timeloop")
    return nothing
end
simulation(data="data/$(Dates.today())curved_drop_test_00")
