# Tutorials

*Swalbe.jl* is software that simulates the dynamics of a thin liquid film.
At it's core is the lattice Boltzmann method (LBM) a mesoscale method that works based on particle distribtuion functions $f_i(\mathbf{x},t)$.
To get a better understandig of both the LBM and the package we supply instructive simulations and discuss the content of the code.

## Constant interface

A nice and easy way to check the consistency of a fluid mechanics solver is to prob the mass or density.
Per continuity equation $\partial_t\rho + \nabla\cdot\mathbf{j} = 0$, we require that the local density $\rho$ can only change if there is a flux $\mathbf{j}$ (okay a divergence of a flux). 
In case of vanishing or constant flux $\nabla\cdot\mathbf{j} = 0$ the density has to be time independent and therefore $\partial_t \rho = 0$. 

Thin film flows can be described with the same kind of equation

```math
\partial_t h(\mathbf{x},t) + \nabla\left(M(h)\nabla p(\mathbf{x},t)\right) = 0. 
```

The local height can only change with time if some pressure gradient 

Having a vanshing pressure gradient the height $h(\mathbf{x},t)$ have to be constant independent of time.
To check this we put together a small sample simulation that probes this.

```julia
using Swalbe # import the package

# Define the system size and parameters
sys = Swalbe.SysConst(Lx = 100,     # 100 lattice units in x-direction 
                      Ly = 100,     # 100 lattice units in y-direction
                      Tmax = 1000)  # LBM time loop runs for 1000 iterations


# Allocation of distribution functions, macroscopic variables and forces 
fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys,      # Size 
                                                                                                          "CPU",    # Where to run 
                                                                                                          false,    # No fluctuations  
                                                                                                          Float64)  # Num type
# Constant film thickness, or flat interface
height .= 1.0
# Compute the first equilibrium distribtuion function based on the initial conditions 
Swalbe.equilibrium!(ftemp, height, velx, vely, vsq)
# Empty list to measure the total mass over time
mass = []
# To add the forcing due to a pressure gradient set true
pressure_gradient = false
# Start of the lattice Boltzmann time loop
for t in 1:sys.Tmax
    # Fill the list 
    push!(mass,sum(height))
    # Talks with us all t % tdump time sets
    if t % sys.tdump == 0
        println("Time step $t mass is $(round(sum(height), digits=3))")
    end
    if pressure_gradient
        # The pressure computation 
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        # Pressure gradient
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        # Slippage
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        # Summation of forces
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
    end
    # New equilibrium
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    # Collision and streaming of distribution function
    Swalbe.BGKandStream!(fout, feq, ftemp, Fx, Fy)
    # Update of macroscopic quantities
    Swalbe.moments!(height, velx, vely, fout)
    # Next iteration
end
# Check if the mass is constant
using Plots
# Since h = 1 everywhere and the size 100 x 100 the mass has to be 100^2
plot(mass, xlabel="time [Δt]", ylabel="Mass [lbu]", label="flat interface", ylim=(100^2-0.1, 100^2+0.1))
```

In fact this is a particular strength of the lattice Boltzmann method.
Under the assumption that no forces are applied the mass is mathematically conserved.
Which is shown in the lower plot

![flat](https://user-images.githubusercontent.com/26249811/124790156-39313700-df4b-11eb-93af-9957bb46e411.png)

## Dewetting of patterned a substrate

Next we take a look at substrates with a wettability gradient and show how to use an image or geometrical shapes to induce directed dewetting.
Here we actually reuse the display simulation of the [README.md](https://github.com/Zitzeronion/Swalbe.jl/blob/master/README.md).
The idea is to use [Youngs equation](https://upload.wikimedia.org/wikipedia/commons/8/85/Thomas_Young-An_Essay_on_the_Cohesion_of_Fluids.pdf)

```math
\cos(\theta_{\text{eq}}) = \frac{\gamma_{sg} - \gamma_{sl}}{\gamma}, 
```

where the $\gamma$'s are the three interfacial energies and $\theta_{\text{eq}}$ is the equilibrium contact angle to address the wettability of the substrate.
Without problems it is possible to discretize $\theta_{\text{eq}}$ similar to the pressure or the film thickness and therefore effectively introduce a wettability gradient $\nabla \theta_{\text{eq}} \neq 0$.
In the logo simulation what actually happens is that initial perturbations `(ϵ * randn())` grow with time, leading to a film rupture and travelling rims.
The ruptures are triggered in regions of high contact angle and the rims meet in regions of low contact angle (letters of the logo).
There is a lot of theory involved and if you want to read up on it check out the review from [Craster&Matar](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131) (paywall), especially section 5. Films driven by intermolecular forces.

But let's start again with some code, first with a triangle with lower wettability in the middle
```julia
using Swalbe

# Define the system size and parameters
sys = Swalbe.SysConst(Lx = 100,     # 100 lattice units in x-direction 
                      Ly = 100,     # 100 lattice units in y-direction
                      n = 3,        # first disjoining pressure exponent
                      m = 2,        # second disjoining pressure exponent
                      Tmax = 1000)  # LBM time loop runs for 1000 iterations


# Allocation of distribution functions, macroscopic variables and forces 
fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys,      # Size 
                                                                                                          "CPU",    # Where to run 
                                                                                                          false,    # No fluctuations
                                                                                                          Float64)  # Num type
# Random (gaussian) perturbations put on top of a flat interface
ϵ = 0.001
height .= 1.0 .+ ϵ .* randn(sys.Lx, sys.Ly)
# The contact angle field
theta = fill(1/6, sys.Lx, sys.Ly)
# Lower contact angle in the middle of the substrate
Swalbe.trianglepattern(theta, 1/6, δₐ=-1/18)
# Compute the first equilibrium distribution function based on the initial conditions 
Swalbe.equilibrium!(ftemp, height, velx, vely, vsq)
# Difference in total height with time
diff_h = []
# Start of the lattice Boltzmann time loop
for t in 1:sys.Tmax
    # Fill the difference list
    push!(diff_h, maximum(height) - minimum(height))
    # Talks with us all t % tdump time sets
    if t % sys.tdump == 0
        println("Time step $t mass is $(round(sum(height), digits=3))")
    end
    # The pressure computation 
    Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, theta, sys.n, sys.m, sys.hmin, sys.hcrit)
    # Pressure gradient
    Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
    # Slippage
    Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
    # Summation of forces
    Fx .= -h∇px .- slipx
    Fy .= -h∇py .- slipy
    
    # New equilibrium
    Swalbe.equilibrium!(feq, height, velx, vely, vsq)
    # Collision and streaming of distribution function
    Swalbe.BGKandStream!(fout, feq, ftemp, Fx, Fy)
    # Update of macroscopic quantities
    Swalbe.moments!(height, velx, vely, fout)
    # Next iteration
end
# Another library for plotting and to my understanding actually the best you can do
using CairoMakie
# We are interested in the heatmap of the film thickness and the growth rate of the perturbation
let
	x1 = 1:1:sys.Tmax
	y2 = zeros(sys.Tmax)
	y2 .= diff_h
    fig = Figure(resolution = (960,450))
	
    ax1 = Axis(fig, xlabel = "time [Δt]", ylabel = "Δh [lbu]",xscale = log10, yscale = log10,xgridstyle=:dash, ygridstyle=:dash, xminorticksvisible = true,
        xminorticks = IntervalsBetween(9), yminorticksvisible = true,
        yminorticks = IntervalsBetween(9))
	leg = lines!(ax1, x1, y2, color = :navy)
	ax2 = Axis(fig,  aspect = 1, xlabel = "x [Δx]", ylabel = "y [Δx]")
    hmap = heatmap!(ax2, height, colormap = :viridis)
    cbar = Colorbar(fig, hmap, label = "thickness", ticksize=15, tickalign = 1, width = 15)
	fig[1,1] = ax1
    fig[1,2] = ax2
    fig[1,3] = cbar
    fig
end
```

Fluid is drained into regions of lower contact angle, therefore into in the triangle.
The effect is the strongest around the vertices of the triangle.
Since in principle this a dewetting instability with a *well defined spectrum* computable using the surface tension $\gamma$ and the disjoining pressure functional $\Pi(h)$.
The result should look like the plot below 

![patterned](https://user-images.githubusercontent.com/26249811/125091450-b3d78f00-e0d0-11eb-9786-9a817a495139.png)

Of course you can and should play with the parameters ($\gamma$, $\delta$, $h_0$, ...) to get a physically correct simulation :wink:.

## Droplet spreading in 1D

If a certain amount of liquid is placed on a surface, like a rain drop :droplet: hitting a plant leaf :four_leaf_clover: we observe a drop sticking to the leaf.
The shape of the droplet is actually nature solving Young's law and finds an equilibrium shape which can be described by a simple macroscopic observable, the contact angle $\theta_{\text{eq}}$.
There is of course more to it, for further reading check out [Snoeijer & Andreotti](http://stilton.tnw.utwente.nl/people/snoeijer/Papers/2013/SnoeijerARFM13.pdf)
We can recast this behavior with a simple and fast simulation.

In contrast to the two examples before we put the whole experiment into a function and just call the function.
While this may seem overkill here it is in fact very useful.
Instead of writing a script for every experiment, we simply write a function and loop through function arguments.
Making it very convenient to perform phase space scans and parameter studies.

```julia
using Swalbe
# Simulation that let's a droplet relax towards it's equilibrium contact angle
function run_dropletrelax(
    sys::Swalbe.SysConst_1D;    # System Constants
    radius=20,                  # Initial droplet radius
    θ₀=1/6,                     # Initial contact angle
    center=(sys.L÷2),           # Position of the center of mass
    verbose=true,               # Simulation prints output 
    T=Float64                   # Number subtype
)
    println("Simulating an out of equilibrium droplet")
    # Empty list to store the radius evolution
    diameter = []
    # Allocation
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p = Swalbe.Sys(sys, false, T)
    # Initial condition, see initialvalues.jl, or ?Swalbe.singledroplet
    Swalbe.singledroplet(height, radius, θ₀, center)
    # Initial equilibrium, in this case a D1Q3 equilibrium
    Swalbe.equilibrium!(ftemp, height, vel)
    # Lattice Boltzmann loop starts
    for t in 1:sys.Tmax
        if verbose
        # Check if the mass conserved
            if t % sys.tdump == 0
                mass = 0.0
                mass = sum(height)

                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        # Push the number of lattice sides inside the droplet to the list
        push!(diameter, length(findall(height .> 0.06)))
        # Compute film pressure with contact angle \theta = 1/9
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
        # Compute the gradient of the pressure and multiply it with the height
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        # Calculate the substrate friction, velocity boundary condition
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        # Sum the forces up
        F .= h∇p .+ slip
        # New equilibrium
        Swalbe.equilibrium!(feq, height, vel)
        # Collide and stream
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        # Compute the new moments
        Swalbe.moments!(height, vel, fout)
    end
    return height, diameter
end
```

If we defined the function as above we can use it to run several experiments and test for example [Tanners law](https://iopscience.iop.org/article/10.1088/0022-3727/12/9/009/pdf).
The idea of Tanner was that the evolution of the droplets radius during spreading should be captured by a powerlaw $R(t) \propto t^{\alpha}$, with $\alpha = 1/10$ in this case.
There are subtleties to this which you can read up in this very nice paper by [Eddi et al.](http://stilton.tnw.utwente.nl/people/eddi/Papers/PhysFluids_spreading.pdf).
**Swalbe.jl** by definition of a numerical solver does not know about real world experiments.
That is why we have to find the correct parameters to capture experimental findings (real world physics), in this case we like to observe a powerlaw growth in radius with $\alpha = 1/10$. 
There are two things we could easily change, the surface tension $\gamma$ and the velocity boundary or slippage $\delta$.

```julia
# Dictionary to store the results
results = Dict()
# Loop over different slip lengths
for slip in [2.0, 1.0, 0.5]
    # Simulation parameters
    sys = Swalbe.SysConst_1D(L=2048, γ=0.001, n=3, m=2, δ=slip, Tmax=2000000);
    # Run the simulation
    h, d = run_dropletrelax(sys, radius=400, θ₀=1/4)
    # Store the data in the dict
    results[slip] = d
end
```

Given the finite slippage we do not observe *large* deviations from the $\alpha=1/10$ powerlaw in the long time limit.


## Active film

In [Richter et all](https://arxiv.org/abs/2402.14635) we introduced a thin film model for chemical reactions inside a liquid film aimed at Supported liquid phase catalysts. For the details of the model look at the article. Here we give a minimal example of how to run a simulation, starting from a homogenous film, leading to film rupture. 

Lets get straight in it 
```julia 
using Swalbe, DelimitedFiles
function simulation(;
    D=5e-6,			# Diffusion coefficient
    Tmax=1e3, 			# Simulation time
    tdump=Int(floor(Tmax/100)), # dump time
    h= 1,  			# initial film height
    rho=1,			# initial catalyst concentration
    ϵ=0.001,			# initial noise
    γ_0=0.001,			# reference surface tension in the absence of chemical compounds
    Γ = 0.1, 			# Surface tension effect product rho_A
    GammaB=0.1,			# Surface tension effect reactant rhoB
    data = "data",		# location to save data
    L=2^12, 			# system size
    alpha=0.0, 			# catalyst vertical distribution parameter 
    theta_0=1/9,		# reference contact angle in fractions of pi
    n=9, m=3, hmin=0.1, hcrit=0.05, hcrit2=0.05,	# disjoining pressure parameters
    delta=1.0,			# slip length
    production_rate=0.1,	# production rate
    sigma_A_up=0.1, 		# evaporation rate product rho_A
    sigma_B_up=0.1, 		# evaporation rate reactant rho_B
    rho_BRes_sigma_B_down=0.1,	# sorption rate reactant rho_B
    rho_0B=rho_BRes_sigma_B_down*h/(production_rate*rho + sigma_B_up),	# initial reactant concentration 
    rho_0A=production_rate*rho*rho_0B/sigma_A_up,			# initial product cocentration 
    D_Product=D*20, D_B=D*20,	# Diffusion coefficients Product rho_A and reactant rho_B
    rho_crit=0.05,		# precursor length catalyst
    mu=1/6,			# viscosity
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
        Swalbe.browniannoise!(state.height, h, ϵ, L/2-L/10)
	i = 0
	# The LBM time loop
	for t in 0:sys.Tmax
		# dump data
		if t % sys.tdump == 0 
			writedlm("$(data)/$(Swalbe.to_4_digits(i))_height.csv", state.height)
			writedlm("$(data)/$(Swalbe.to_4_digits(i))_rho.csv", state.rho)
			writedlm("$(data)/$(Swalbe.to_4_digits(i))_rho_A.csv", state.rho_A)
			writedlm("$(data)/$(Swalbe.to_4_digits(i))_rho_B.csv", state.rho_B)
			i += 1
			println("Time step t=$t")
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
simulation()
```


That's it

## Thin films on curved substrate

It is possible to include a curved substrate into the model using the substrate profile `b` into the pressure
```math
p_{\text{cap}} = -(1 + 0.5 \nabla s^2) \gamma (1 - \cos(\theta)) \phi'(h) - \gamma \nabla^2 (h+b)
```

Here is an example on how to do this 

```julia 
using Swalbe, DelimitedFiles, Dates
function simulation(;
	Tmax=1e8, 			# Simulation time
	tdump=Int(floor(Tmax/100)),	# dump time
	h= 1,				# initial film height  
	eps=0.001, 			# initial noise
	gamma=1/6^2,			# surface tension
	data = "data",			# folder for saving data
 	L=2^12,				# system size
    	theta=1/9,			# contact angle
    	n=9, m=3, hmin=0.1, 		# parameters for disjoining pressure
	delta=1.0,			# slip length, in combination with heigher hmin i
	omega=3,			# number of periods of the underlying substrate
)
	# set up system
    	sys = Swalbe.SysConst_1D{Float64}(L =  L ,Tmax=Tmax, tdump=tdump, γ = gamma, δ=delta,θ =theta, n=n,m=m,hmin=hmin)
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
    	for t in 0:sys.Tmax
		# save data
        	if t % sys.tdump == 0
                	writedlm("$(data)/$(Swalbe.to_4_digits(i))_height.csv", state.height)
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
simulation()
```


## Miscible liquids

Here I show to simulate the colaescence of two droplets of miscible liquids with different surface tensions, lets say water and ethanol. 

```julia
using Swalbe, DelimitedFiles,  Dates


function run_Miscible(;
    eps=0.01,                                   # initial noise
    Tmax=Int(1e7),                              # Simulation time
    dumps=100,                                  # number of dumps
    tdump=Int(floor(Tmax/dumps)),               # dump time
    L=2^9,                                      # system size
    delta=1.0,                                  # slip length
    hmin=0.1, hcrit=0.05, n=9, m=3,             # disjoining pressure parameters
    gamma= 0.01 .* ones(3,2),                   # surface tensions
    data="path/to/data",                        # save path
    D=1e-2,                                     # diffusion coeficient
    r=50,                                       # Initial droplet radiusj
    center_weight=0.5                           # weight for surface tension smoothing
    )
    #Make folders
    datathis= "$(data)/Tmax=$Tmax,L=$L,delta=$delta,gamma=$(replace(replace("$gamma", ";"=>","), " " => "_")),D=$D,r=$r,hmin=$hmin,hcrit=$hcrit"
        try
        mkdir(dataall)
        catch
        end
        try
        mkdir(datathis)
        catch
        end
        try
            mkdir("$(datathis)/figs")
        catch
        end
    #Setup system
    sys = Swalbe.SysConstMiscible_1D(L=L, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=delta, n=n, m=m, D=D, hmin=hmin, hcrit=hcrit)
    state = Swalbe.Sys(sys)
    # Initial data
    theta_1=acos((sys.gamma[3,1]-sys.gamma[2,1])/sys.gamma[1,1])
    theta_2=acos((sys.gamma[3,2]-sys.gamma[2,2])/sys.gamma[1,2])
    state.height[:,1] .= Swalbe.droplet_base(state.height[:,1], r, theta_1/pi, Int(sys.L/2 - r/2), precursor=(sys.hmin-2*sys.hcrit)/2)
    state.height[:,2] .= Swalbe.droplet_base(state.height[:,2], r, theta_2/pi, Int(sys.L/2 + r/2), precursor=(sys.hmin-2*sys.hcrit)/2)
    bridge_height=[]
    i=0
    println("Starting LBM time loop")
    for t in 0:sys.Tmax
        if t % sys.tdump == 0
            #Store data
            writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_h_t=$t",state.height)
            println("time step t=$t")
            i+=1
        end
        #actual simulations
            # surface tension calculation from the height fields
            # Swalbe.surface_tension!(state,sys)
            # surface tension calculation from the height fields with moving averadge smoothing applied
            Swalbe.surface_tension_smooth!(state,sys, center_weight=center_weight)
            Swalbe.filmpressure!(state, sys)
            Swalbe.h∇p!(state,sys)
            Swalbe.slippage!(state,sys)
            # calculate Marangoni shear stress
            Swalbe.∇gamma!(state, sys)
            state.F .= -state.h∇p  .-  state.slip .+ state.stress
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state, sys)
            Swalbe.moments!(state, sys)
    end
end


dataall = "data/$(Dates.today())Miscible_coal_test_00"
gamma=ones(3,2)
gamma_0=0.002
# fixed contact angles
Delta=0.2
gamma[1,1]=gamma_0*(1+Delta/2)
gamma[1,2]=gamma_0*(1-Delta/2)
# solid vapour
gamma[3,1]=2*gamma_0
gamma[3,2]=2*gamma_0
theta=1/9
# fixed contact angle
gamma[2,1]=gamma[3,1]-gamma[1,1]*cospi(theta)
gamma[2,2]=gamma[3,2]-gamma[1,2]*cospi(theta)
run_Miscible(Tmax=Int(1e7),gamma=gamma, D=1.5e-4, r=200, delta=1.0, hmin=0.2, L=2^10, data=dataall)
```

Have fun with it

## Multilayer thin films

This shows how to simulate a layer film resting on a layer of water. This we call a two layer multilayer system. To read more about the theory see [Richter et al.](https://arxiv.org/abs/2409.16659). 

This is the code on how to do this: 
```julia 
using Pkg

Pkg.activate(".")


using Swalbe, DelimitedFiles, CUDA


function run_multilayer(;
    h_1=2,		# initial height lower layer
    h_2=1,		# initial height upper layer
    eps=[0.001,0.001],	# initial noise
    Tmax=Int(1e7),	# simulation time
    dumps=100,		# number of dumps
    tdump=Int(floor(Tmax/dumps)), 	# dump time
    Lx=2^8, Ly=2^8,	# system size
    delta=1.0,		# slip lower layer
    delta_2=1.0,	# slip upper layer
    n=9, m=3,		# disjoining pressure exponents
    gamma= [0 0.001 0.001 0.001 ; 0 0 0.001 0.001 ; 0 0 0 0.001; 0 0 0 0],	# surface tensions gamma[i,j] is the interaction between layer i-1 and j-1 having 0 for the solid substrat and 3 the atmosphere
    data="data", 	# path to save simulation dumps
    hmin=0.1,		# precursor layer
    mu=[1/6 1/6],	# viscosities 
    )
    #Setup system
    sys = Swalbe.SysConstMultiLayer(Lx=Lx, Ly=Ly, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=[delta,delta_2], n=n, m=m, hmin=hmin, mu=mu)
    state = Swalbe.Sys(sys, device=device)
    # this is needed to send back and force in between GPU and CPU
    host_height=zeros(Lx,Ly,2)
    #initial conditions
    host_height[:,:,1] .= Swalbe.browniannoise2D!(host_height[:,:,1], h_1, eps[1], Lx/2-Lx/10)
    host_height[:,:,2] .= Swalbe.browniannoise2D!(host_height[:,:,2], h_2, eps[2], Lx/2-Lx/10)
    CUDA.copyto!(state.height, host_height)
    i=0
    println("Starting LBM time loop")
    for t in 0:sys.Tmax
        if t % sys.tdump == 0 
            #Store data
                CUDA.copyto!(host_height,state.height)
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_h1_t=$t",host_height[:,:,1])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_h2_t=$t",host_height[:,:,2])
            #status report
            println("Time step t=$t)
            i+=1
        end
        #actual simulations
        if run 
            Swalbe.filmpressure!(state, sys)
	    Swalbe.h∇p!(state,sys)
            Swalbe.slippage!(state, sys)
	    state.Fx .= -state.h∇px  .+  state.slipx
	    state.Fy .= -state.h∇py  .+  state.slipy
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state)
            Swalbe.moments!(state, sys)
        end
    end 
end


dataall = "path\to\save\dumps"
gamma_0=0.01
gamma=gamma_0*ones(4,4)
gamma[2,4]=gamma_0*(1+cospi(2/9))
gamma[1,3]=gamma_0*(2.1)
gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]
run_multilayer(Tmax=Int(1e6),gamma=gamma, h_2=1.0, h_1=6.0, delta=1.0, delta_2=1.0, Lx=2^8,Ly=2^8, mu=[0.1 0.1], hmin=0.2, devicenumber=2)
```

## Further tutorials

More tutorials will follow in the future.
I plan to create one for every  paper the method was used for.
So be sure to check out the docs every now and then.

The next tutorial will be about switchable substrates. 
In this case the wettability can not only addressed locally but also with a time dependency.
Here is what happens if the time frequency is [high](https://gist.github.com/Zitzeronion/116d87978ece82c8ae64a3f7edb9dbb3#gistcomment-3811414) and this happens if we update with a [lower](https://gist.github.com/Zitzeronion/116d87978ece82c8ae64a3f7edb9dbb3#gistcomment-3811414) frequency.
