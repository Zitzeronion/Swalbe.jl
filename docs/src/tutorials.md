# Tutorials

*Swalbe.jl* is software that simulates the dynamics of a thin liquid film.
At it's core is the lattice Boltzmann method (LBM) a mesoscale method that works based on particle distribtuion functions $f_i(\mathbf{x},t)$.
To get a better understandig of both the LBM and the package we supply instructive simulations and discuss the content of the code.

## Constant interface

A nice and easy way to check the consistency of a fluid mechanics solver is to prob the mass or density.
Per continuity equation $\partial_t\rho + \nabla\cdot\mathbf{j} = 0$, we require that density $\rho$ can only change if there is a flux $\mathbf{j}$ (okay a divergence of a flux). 
In case of vanishing or constant flux $\nabla\cdot\mathbf{j} = 0$ the density has to be time independent and therefore $\partial_t \rho = 0$. 

Thin film flows can be described with the same kind of equation

$$ \partial_t h(\mathbf{x},t) + \nabla\left(M(h)\nabla p(\mathbf{x},t)\right) = 0. $$

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
# Start of the lattice Boltzmann time loop
for t in 1:sys.Tmax
    # Fill the list 
    push!(mass,sum(height))
    # Talks with us all t % tdump time sets
    if t % sys.tdump == 0
        println("Time step $t mass is $(round(sum(height), digits=3))")
    end
    # The pressure computation 
    Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, 1/9, sys.n, sys.m, sys.hmin, sys.hcrit)
    # Pressure gradient
    Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
    # Slippage
    Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
    # Summation of forces
    Fx .= h∇px .+ slipx
    Fy .= h∇py .+ slipy
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
And this is how the result should look like

![flat](https://user-images.githubusercontent.com/26249811/124790156-39313700-df4b-11eb-93af-9957bb46e411.png)

## Dewetting of patterned a substrate

Next we take a look at patterend substrates and show how to use an image or some geometrical shape to induce directed dewetting.
