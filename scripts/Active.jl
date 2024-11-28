using Pkg; Pkg.activate(".")
using Swalbe, DelimitedFiles, Dates, LaTeXStrings
using Plots; Plots.gr(); ENV["GKSwstype"] = "100"
function simulation(;
    D=5e-6,                     # Diffusion coefficient
    Tmax=1e7,                   # Simulation time
    tdump=Int(floor(Tmax/100)), # dump time
    h= 1,                       # initial film height
    rho=1,                      # initial catalyst concentration
    ϵ=0.001,                    # initial noise
    γ_0=0.001,                  # reference surface tension in the absence of chemical compounds
    Γ = 0.1,                    # Surface tension effect product rho_A
    GammaB=0.1,                 # Surface tension effect reactant rhoB
    data = "data",              # location to save data
    L=2^12,                     # system size
    alpha=0.0,                  # catalyst vertical distribution parameter
    theta_0=1/9,            	# reference contact angle
    n=9, m=3, hmin=0.1, hcrit=0.05, hcrit2=0.05,        # disjoining pressure parameters
    delta=1.0,                  # slip length
    production_rate=0.1,        # production rate
    sigma_A_up=0.1,             # evaporation rate product rho_A
    sigma_B_up=0.1,             # evaporation rate reactant rho_B
    rho_BRes_sigma_B_down=0.1,  # sorption rate reactant rho_B
    rho_0B=rho_BRes_sigma_B_down*h/(production_rate*rho + sigma_B_up),  # initial reactant concentration
    rho_0A=production_rate*rho*rho_0B/sigma_A_up,                       # initial product cocentration
    D_Product=D*20, D_B=D*20,   # Diffusion coefficients Product rho_A and reactant rho_B
    rho_crit=0.05,              # precursor length catalyst
    mu=1/6,                     # viscosity
    plot=true			# do plots
)
	    # optionally create the data folder
	    try 
		mkdir(data)
	    catch e println(e) end
            # Set up sys consts
    	    GammafacA=rho_0A/(h*γ_0)
    	    GammafacB=rho_0B/(h*γ_0)
            sys = Swalbe.SysConstActive_1D{Float64}(
                L =  L,
                D=D,
                Γ = Γ/GammafacA,
                GammaB=GammaB/GammafacB,
                Tmax=Tmax,
                tdump=tdump,
                γ_0 = γ_0,
                δ=delta,
                alpha=alpha,
                θ_0 =theta_0,
                n=n,m=m,hmin=hmin,hcrit=hcrit, hcrit2=hcrit2,
                production_rate=production_rate,
                sigma_A_up=sigma_A_up,
                sigma_B_up=sigma_B_up,
                rho_BRes_sigma_B_down=rho_BRes_sigma_B_down,
                D_Product=D_Product,
                D_B=D_B,
                rho_crit=rho_crit,
                b_A=rho_crit,
                b_B=rho_crit,
                μ=mu
        )
        # set up simulation fields
        state = Swalbe.Sys(sys)
        # Set initial values
        state.rho.=rho
        state.rho_A .= rho_0A
        state.rho_B .= rho_0B
        # This set up the film height as constant plus a small noise without high wavenumber frequencies
        Swalbe.browniannoise!(state.height, h, ϵ, L/2-L/10)
        # The LBM time loop
	i=0
        for t in 0:sys.Tmax
                # dump data
                if t % sys.tdump == 0
                        writedlm("$(data)/$(Swalbe.to_4_digits(i))_height.csv", state.height)
                        writedlm("$(data)/$(Swalbe.to_4_digits(i))_rho.csv", state.rho)
                        writedlm("$(data)/$(Swalbe.to_4_digits(i))_rho_A.csv", state.rho_A)
                        writedlm("$(data)/$(Swalbe.to_4_digits(i))_rho_B.csv", state.rho_B)
			println("time step t=$t, mass=$(sum(state.height)), diff=$(maximum(state.height)-minimum(state.height))")
			if plot
				Plots.plot(state.height, label=L"h")
				Plots.plot!(state.rho, label=L"\rho", linestyle=:dash)
				Plots.plot!(state.rho_A, label=L"\rho_{Product}", linestyle=:dash)
				Plots.plot!(state.rho_B, label=L"\rho_{Reactant}", linestyle=:dash)
			end
			i += 1
                end
                # system update, mostly as in a standard TFE simulation, the parts that are different are commented
                Swalbe.surface_tension_4_fields!(state, sys)
                Swalbe.filmpressure_fast!(state, sys)
                Swalbe.h∇p!(state)
                Swalbe.slippage!(state, sys)
                # claculate the Marangoni stress
                Swalbe.∇γ!(state, sys)
                # Include the Marangoni stress into the forcing
                state.F .= -state.h∇p .- state.slip .+ state.stress
                # Do the LBM loops for catalys state.rho, product state.rho_A, and reactant state.rho_B
                Swalbe.update_rho_LBM!(state,sys)
                Swalbe.update_rho_A_LBM_precursor!(state,sys)
                Swalbe.update_rho_B_LBM_precursor!(state,sys)
                Swalbe.equilibrium!(state)
                Swalbe.BGKandStream!(state, sys)
                Swalbe.moments!(state)
        end
end

# run the simulation
simulation(data="data/$(Dates.today())Active_test/")
