# Tutorials

*Swalbe.jl* is software that simulates the dynamics of a thin liquid film.
At it's core is the lattice Boltzmann method (LBM) a mesoscale method that works based on particle distribtuion functions $f_i(\mathbf{x},t)$.
To get a better understandig of both the LBM and the package we supply instructive simulations and discuss the content of the code.

## Constant interface

A nice and easy way to check the consistency of a fluid mechanics solver is to prob the mass or density.
Per continuity equation $\partial_t\rho + \nabla\cdot\mathbf{j} = 0$, we require that the local density $\rho$ can only change if there is a flux $\mathbf{j}$ (okay a divergence of a flux). 
In case of vanishing or constant flux $\nabla\cdot\mathbf{j} = 0$ the density has to be time independent and therefore $\partial_t \rho = 0$. 

Thin film flows can be described with the same kind of equation

$$ \partial_t h(\mathbf{x},t) + \nabla\left(M(h)\nabla p(\mathbf{x},t)\right) = 0. $$

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
                                                                                                          false,    #      
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

$$ \cos(\theta_{\text{eq}}) = \frac{\gamma_{sg} - \gamma_{sl}}{\gamma}, $$

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
                                                                                                          false,    #      
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
p1 = plot(1:sys.Tmax+1, diff_h, axis=:log, xlabel="time [Δt]", ylabel="thickness difference [lbu]", label="Δh")
p2 = heatmap(height, aspect_ratio=1, xlabel="x [Δx]", ylabel="y [Δx]", c=:viridis, colorbar_title = "film thickness [lbu]" )

plot(p1,p2)
```
