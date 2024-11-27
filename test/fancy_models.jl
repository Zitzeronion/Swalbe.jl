# This tests merely if all model run error free from a numerical point of view. This does not include physical tests or unit tests of the underlying models.
# (C) Tilman
@testset "Fancy Models" begin 
	@testset "Catalytic Thin Film" begin 
	println("Testing catalytic thin films")
	function catalysis(;
    		D=5e-6,                     # Diffusion coefficient
    		Tmax=1e3,                   # Simulation time
    		tdump=Int(floor(Tmax/10)), # dump time
    		h= 1,                       # initial film height
    		rho=1,                      # initial catalyst concentration
    		eps=0.001,                    # initial noise
    		γ_0=0.001,                  # reference surface tension in the absence of chemical compounds
    		Γ = 0.1,                    # Surface tension effect product rho_A
    		GammaB=0.1,                 # Surface tension effect reactant rhoB
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
		)
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
			state.height .= h .* ones(L) .* eps .* rand(L)
			height_0 = zeros(L)
			height_0 .= state.height
        			# The LBM time loop
        		i=0
        		for t in 0:sys.Tmax
				if t%sys.tdump ==0
					println("timestep t=$t")
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
			return height_0, state.height
		end
		height_0, height_end = catalysis()
		@test isapprox(sum(height_0), sum(height_end))
	end
	@testset "Curved Substrates" begin 
	println("testing liquid films on curved substrates")
		function curved(;
        		Tmax=1000,                       # Simulation time
        		tdump=Int(floor(Tmax/10)),     # dump time
        		h= 1,                           # initial film height
        		eps=0.001,                      # initial noise
        		gamma=1/6^2,                    # surface tension
        		L=2^12,                         # system size
        		theta=1/9,                      # contact angle
        		n=9, m=3, hmin=0.1,             # parameters for disjoining pressure
        		delta=1.0,                      # slip length, in combination with heigher hmin i
        		omega=3,                        # number of periods of the underlying substrate
		)
        		# set up system
			sys=Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(Tmax=Tmax, tdump=tdump, γ = gamma, δ=delta,θ =theta, n=n,m=m,hmin=hmin))
        		state = Swalbe.Sys(sys, kind="curved")
        		# Initial data
        		state.height.= h .+ eps .* rand(sys.L)
			height_0=zeros(L)
			height_0.=state.height
        		# set up system profile
        		state.substrate .=[1/4*sin(i*omega*2*pi/sys.L) for i in 1:sys.L]
        		# We need the gradient of the substrate profile
        		Swalbe.∇f!(state.grad_substrate, state.substrate, state.dgrad)
        		println("Starting Lattice Boltzmann time loop for active thin film")
        		i=0
        		for t in 0:sys.Tmax
				if t%sys.tdump == 0 
					println("timestep t=$t")
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
    		return height_0, state.height
		end
		height_0, height_end = curved()
		@test isapprox(sum(height_0), sum(height_end))
	end
	@testset "Miscible liquid films" begin 
	println("Testing miscible liquid films")
		function run_Miscible(;
    			eps=0.01,
    			Tmax=Int(1e3),
    			dumps=10,
    			tdump=Int(floor(Tmax/dumps)),
    			L=2^9,
    			delta=1.0,
    			hmin=0.1,
    			hcrit=0.05,
    			n=9, m=3,
    			gamma= [0.002 0.002; 0.00202702 0.00202702; 0.004 0.004],
    			tau=1.0,
    			D=1e-2,
    			h_2=3,
    			r_1=50,
    			center_weight=0.5
    			)
    			sys = Swalbe.SysConstMiscible_1D(L=L, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=delta, n=n, m=m, tau=tau, D=D, hmin=hmin, hcrit=hcrit)
    			state = Swalbe.Sys(sys)
    			theta_1=acos((sys.gamma[3,1]-sys.gamma[2,1])/sys.gamma[1,1])
    			theta_2=acos((sys.gamma[3,2]-sys.gamma[2,2])/sys.gamma[1,2])
    			bar_theta=(theta_1+theta_2)/2
    			println("theta_1=$theta_1, theta_2=$theta_2, bar_theta=$bar_theta")
    			delta_gamma=abs(sys.gamma[1,1]-sys.gamma[1,2])
    			bar_gamma=(sys.gamma[1,1]+sys.gamma[1,2])/2
    			println("Delta_gamma=$delta_gamma, bar_gamma=$bar_gamma")
    			M=3*delta_gamma/(2*bar_gamma*bar_theta^2)
    			println("M=$M")
    			state.height[:,1] .= Swalbe.droplet_base(state.height[:,1], r_1, theta_1/pi, Int(sys.L/2 - r_1/2), precursor=(sys.hmin-2*sys.hcrit)/2)
    			state.height[:,2] .= Swalbe.droplet_base(state.height[:,2], r_1, theta_2/pi, Int(sys.L/2 + r_1/2), precursor=(sys.hmin-2*sys.hcrit)/2)
			height_0=zeros(L,2)
			height_0 .= state.height
    			println("Starting LBM time loop")
    			for t in 0:sys.Tmax
        			if t % sys.tdump == 0
            				println("Time step t=$t, mass=$(sum(state.height))")
				end
				Swalbe.surface_tension_smooth!(state,sys, center_weight=center_weight)
            			Swalbe.filmpressure!(state, sys)
            			Swalbe.h∇p!(state,sys)
            			Swalbe.slippage!(state,sys)
            			Swalbe.∇gamma!(state, sys)
            			state.F .= -state.h∇p  .-  state.slip .+ state.stress
            			Swalbe.equilibrium!(state)
            			Swalbe.BGKandStream!(state, sys)
            			Swalbe.moments!(state, sys)
			end
			return height_0, state.height
		end
		height_0, height_end = run_Miscible(D=1.5e-4, h_2=5, r_1=200, delta=1.0, hmin=0.2, L=2^10)
		@test isapprox(sum(height_0), sum(height_end))
	end
	@testset "2 layer thin film 2D" begin
	println("Testing 2 Layer Thin Films in 2D")
		function run_multilayer(;
    			h_1=2,
    			h_2=1,
    			eps=[0.001,0.001],
    			Tmax=Int(1e3),
    			dumps=10,
    			tdump=Int(floor(Tmax/dumps)),
    			L=2^8,
    			delta=1.0,
    			delta_2=1.0,
    			n=9, m=3,
    			gamma= [0 0.001 0.001 0.001 ; 0 0 0.001 0.001 ; 0 0 0 0.001; 0 0 0 0],
    			plot=true,
    			run=true,
    			hmin=0.1,
    			hcrit=0.05,
    			mu=[1/6 1/6]
    			)
			sys = Swalbe.SysConstMultiLayer_1D(L=L, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=[delta,delta_2], n=n, m=m, hmin=hmin, hcrit=hcrit, mu=mu)
    			state = Swalbe.Sys(sys)
 			state.height[:,1] .= h_1 .* ones(L) .+ eps[1] .* rand(L)
 			state.height[:,2] .= h_2 .* ones(L) .+ eps[2] .* rand(L)
			height_0 = zeros(L,2)
			height_0 .= state.height
			for t in 0:sys.Tmax
        			if t % sys.tdump == 0
					println("time step t=$t")
				end
            			Swalbe.filmpressure!(state, sys)
            			Swalbe.h∇p!(state,sys)
				Swalbe.slippage!(state,sys)
                		state.F .= -state.h∇p  .+  state.slip
				Swalbe.equilibrium!(state)
            			Swalbe.BGKandStream!(state)
            			Swalbe.moments!(state, sys)
			end
			return height_0, state.height
		end
		height_0, height_end = run_multilayer()
		@test isapprox(sum(height_0), sum(height_end))
	end
	@testset "3 layer thin film 2D" begin 
		println("Testing 3 Layer Thin Films in 2D")
		function run_multilayer3(;
    			h_1=2,
    			h_2=1,
    			h_3=1,
    			eps=[0.001,0.001, 0.001],
    			Tmax=Int(1e3),
    			dumps=10,
    			tdump=Int(floor(Tmax/dumps)),
    			L=2^8,
    			# Only the first entry is actually used
    			delta=[0.0,0.0,0.0],
    			n=9, m=3,
    			gamma= 0.001 .* ones(5,5),
    			plot=true,
    			run=true,
    			hmin=0.1,
    			hcrit=0.05,
    			mu=[1/6 1/6 1/6],
    			extra_pressure_fac=1.0
    			)
    			#Setup system
    			sys = Swalbe.SysConstMultiLayer_1D(L=L, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=delta, n=n, m=m, hmin=hmin, hcrit=hcrit,mu=mu, layers=3, extra_pressure_fac=extra_pressure_fac)
    			state = Swalbe.Sys(sys)
 			state.height[:,1] .= h_1 .+ [0.25*sin(i*2*pi/sys.L) for i in 1:sys.L]
    			state.height[:,2] .= h_2
    			state.height[:,3] .= h_3
			height_0 = zeros(L,3)
			height_0 .= state.height
			for t in 0:sys.Tmax
        			if t % sys.tdump == 0
            				println("Time step t=$t")
				end
            			Swalbe.filmpressure!(state, sys)
            			Swalbe.extra_pressure!(state,sys)
                		Swalbe.h∇p!(state,sys)
           			Swalbe.slippage!(state, sys)
            			# Swalbe.slippage_no_slip!(state, sys)
                		state.F .= -state.h∇p  .+  state.slip
            			Swalbe.equilibrium!(state)
            			Swalbe.BGKandStream!(state)
            			Swalbe.moments!(state, sys)
			end
			return height_0, state.height
        	end
		height_0, height_end = run_multilayer3()
		@test isapprox(sum(height_0), sum(height_end))
	end
	@testset "2 layer thin film in 3D" begin 
		println("testing 2 layers in 3D")
		println("This requires CUDA")
		using CUDA
		function run_multilayer3d(;
    			h_1=2,
    			h_2=1,
    			eps=[0.001,0.001],
    			Tmax=Int(1e3),
    			dumps=10,
    			tdump=Int(floor(Tmax/dumps)),
    			Lx=2^8,
    			Ly=2^8,
    			delta=1.0,
    			delta_2=1.0,
    			n=9, m=3,
    			gamma= [0 0.001 0.001 0.001 ; 0 0 0.001 0.001 ; 0 0 0 0.001; 0 0 0 0],
    			plot=true,
    			run=true,
    			hmin=0.1,
    			mu=[1/6 1/6],
    			device="GPU",
    			devicenumber=1,
    			)
    			sys = Swalbe.SysConstMultiLayer(Lx=Lx, Ly=Ly, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=[delta,delta_2], n=n, m=m, hmin=hmin, mu=mu)
    			state = Swalbe.Sys(sys, device=device)
 			host_height=zeros(Lx,Ly,2)
			host_height[:,:,1] .= h_1 .* ones(Lx, Ly) .+ eps[1] .* rand(Lx,Ly)
			host_height[:,:,2] .= h_2 .* ones(Lx, Ly) .+ eps[2] .* rand(Lx,Ly)
			CUDA.copyto!(state.height, host_height)
    			println("Starting LBM time loop")
    			for t in 0:sys.Tmax
        			if t % sys.tdump == 0
            				println("Time step t=$t")
				end
            			Swalbe.filmpressure!(state, sys)
				Swalbe.h∇p!(state,sys)
            			Swalbe.slippage!(state, sys)
            			state.Fx .= -state.h∇px  .+  state.slipx
            			state.Fy .= -state.h∇py  .+  state.slipy
				Swalbe.equilibrium!(state)
            			Swalbe.BGKandStream!(state)
				Swalbe.moments!(state, sys)
			end
			host_height_end=zeros(Lx,Ly,2)
			CUDA.copyto!(host_height_end, state.height)
			return host_height, host_height_end
		end
		height_0, height_end = run_multilayer3d()
		@test isapprox(sum(height_0), sum(height_end))
	end
end
