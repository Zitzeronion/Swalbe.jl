# Swalbe.jl

A lattice Boltzmann framework to solve thin liquid film problems.

## Thin film simulations using lattice Boltzmann

Why is a thin film solver called **Swalbe.jl** you may ask?

The idea is to use the [*lattice Boltzmann method (LBM)*](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) and all its benefits (easy to code, vast amount of literature and scalability) to simulate thin liquid film flows.
Instead of reinventing the wheel we make use of a class of lattice Boltzmann models that were build to simulate shallow water problems, see [Salmon](http://pordlabs.ucsd.edu/rsalmon/salmon.1999a.pdf) (not the fish :fish:), [Dellar](https://people.maths.ox.ac.uk/dellar/papers/LBshallow.pdf) and [van Thang et al.](https://hal.archives-ouvertes.fr/hal-01625073/document) (*all free to read*).
Thus the name of the package **S**hallow **WA**ter **L**attice **B**oltzmann slov**E**r or **Swalbe**.

Of course using a plain shallow water model will not work to simulate thin film dynamics, that is the reason we build our own model :neckbeard:.
Now the main difference is that we throw away most of the shallow water parts by assuming they are small as compared to thin film relevant things, e.g. the substrate fluid interaction.
The full explanation of the model with some benchmarks can be found in our paper [Zitz et al.](http://pub.hi-ern.de/publications/2019/ZSMDH19/2019-ThinFilm-PRE.pdf) (the C/C++ OpenACC codebase has not been further developed since the project moved to *Julia*)

## How to **get** 

First of all you need a *Julia* installation. 
*Julia* is a high level open source programming language and it is as easy to use as python :snake: (my opinion).

*Julia* can be downloaded at the projects homepage [julialang.org](https://julialang.org/), or clones from the github [repo](https://github.com/JuliaLang/julia). If you download *Julia* from the homepage make sure that you use the correct installation for your operating system.
**Important** for **CUDA** we require *Julia* version 1.6 or higher, usually the most recent version is also the one you should aim for.  

**Swalbe.jl** is a registered package of the *Julia* package manager.
The only thing you have to do is to add the package to your *Julia* environment with: 
```julia
julia> ] add Swalbe
```
Of course you can as well clone or fork the repo and activate the package inside der **REPL**.
First you need to go the Swalbe directory and open a **REPL**
```bash
cd \Swalbe_folder
julia
```

now you can activate the package with
```julia
julia> ] activate .
  Activating environment at `local_Swalbe_folder`
(Swalbe)> 
```

To see that the package works you can run the test suit with
```julia
julia> ] test Swalbe
```
All tests can be found in [test folder](https://github.com/Zitzeronion/Swalbe.jl/tree/master/test), but do not expect too many comments.
Still especially the [simulate.jl](https://github.com/Zitzeronion/Swalbe.jl/blob/master/test/simulate.jl) file is worth a look.

## How to **use** 

The idea of **Swalbe.jl** is to script your thin film simulation, based on a lattice Boltzmann iteration.
That is why most core functions can easily be extended, or used out of the box.
So how does it work, fist we have to load **Swalbe.jl** into the REPL or put the following line on top of our script
```julia
julia> ] add Swalbe         # If not yet added
julia> using Swalbe
```
which can take a minute or so, don't be alarmed if it takes more than ten seconds.
Alternatively one can use [DrWatson](https://github.com/JuliaDynamics/DrWatson.jl) (super cool package to manage scientific computing with *Julia*) and use the `@quickactivate :Swalbe` macro, have a look at the [scripts](https://github.com/Zitzeronion/Swalbe.jl/blob/master/scripts) folder for inspiration.

As a picture says more than a thousand words here is a shiny use case of **Swalbe.jl** 
```julia
using Images, Colors, Swalbe

"""
  dewet_logo(logo_source, kwargs...)

Dewetting of a patterned substrate with pattern according to a image file at `logo_source`.
"""
function dewet_logo(logo_source;         # png file location
                    ϵ=1e-3,              # initial perturbation
                    h₀=1.0,              # initial film thickness
                    device="CPU",        # simulation on the CPU
                    slip=3.0,            # slip length, see three phase contact line
                    Tmax=10000,          # number of lattice Boltzmann iterations
                    dump=100,            # saving interval   
                    T=Float64,           # numeric accuracy
                    verbose=true)        # let's talk         
	println("Starting logo dewetting")
	# Load the image file
    logo = load(logo_source)
    # Set up of the simulation constants
    sys = Swalbe.SysConst(Lx=size(logo)[1], Ly=size(logo)[2], Tmax=Tmax, tdump=dump, δ=slip)
    # Memory allocation
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    # Output
    fluid = zeros(Tmax÷dump, sys.Lx*sys.Ly) 
    theta = zeros(sys.Lx, sys.Ly)
    println("Reading logo: $(logo_source)\nand pattern substrate according to it")
    # Lower contact angle inside the letters, here the red channel of the image is used
    theta = T.(2/9 .- 1/18 .* red.(reverse(rot180(logo), dims=2)))
    if device == "CPU"
        for i in 1:sys.Lx, j in 1:sys.Ly
            # Initial height configuration
            height[i,j] = h₀ + ϵ * randn()
        end
        th = zeros(size(height))
        th .= theta 
    elseif device == "GPU"
        h = zeros(size(height))
        for i in 1:sys.Lx, j in 1:sys.Ly
            h[i,j] = h₀ + ϵ * randn()
        end
        # Lower contact angle inside the letters
		# Forward it to the GPU
		th = CUDA.adapt(CuArray, theta)
        height = CUDA.adapt(CuArray, h)
    end
    # Computation of the initial equilibrium
    Swalbe.equilibrium!(fout, height, velx, vely, vsq)
    ftemp .= fout
    # Lattice Boltzmann time loop
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            deltaH = maximum(height) - minimum(height)
            # Simulation talks with you
            if verbose
                println("Time step $t mass is $(round(mass, digits=3)) and δh is $(round(deltaH, digits=3))")
            end
        end
        # Calculation of the pressure and the pressure gradient
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, th, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        # Forces are the pressure gradient and the slippage due to substrate liquid boundary conditions
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        # New equilibria
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        # Single relaxation and streaming
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        # New moments
        Swalbe.moments!(height, velx, vely, fout)
        # Measurements, in this case only snapshots of simulation's arrays
        Swalbe.snapshot!(fluid, height, t, dumping = dump)
    end
    return fluid
    # Free the GPU
    if device == "GPU"
        CUDA.reclaim()
    end
end
# Run the simulation
dewet_logo("path-to-logo", slip=3.0, Tmax=10000, dump=100, verbose=true)
```

Here the [Images](https://github.com/JuliaImages/Images.jl) and [Colors](https://github.com/JuliaGraphics/Colors.jl) package allow convenient reading of png or jpg files.
To give you an understanding of what happens here we take a look at the different parts.
First of we define our simulation as `function` in this case `dewet_logo()`.
There is one input needed, namely the location of the png or jpg file you want to dewet, in my case I used our institutes [logo](https://gist.github.com/Zitzeronion/807b9a7b2226e65643288df9a8cc1f46/raw/b83988608cd5cdbdb9240e8182050383f442700f/logo_red.png).
Other arguments are keywords that have a default value.

Next step is to define the system we want to simulate, so mostly allocations and initial conditions as well as substrate patterning
```julia
# Load the image file
logo = load(logo_source)
# Set up of the simulation constants
sys = Swalbe.SysConst(Lx=size(logo)[1], Ly=size(logo)[2], Tmax=Tmax, tdump=dump, δ=slip)
# Memory allocation
fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
# Output
fluid = zeros(Tmax÷dump, sys.Lx*sys.Ly) 
theta = zeros(sys.Lx, sys.Ly)
println("Reading logo: $(logo_source)\nand pattern substrate according to it")
# Lower contact angle inside the letters, here the red channel of the image is used
theta = T.(2/9 .- 1/18 .* red.(reverse(rot180(logo), dims=2)))
if device == "CPU"
    for i in 1:sys.Lx, j in 1:sys.Ly
        # Initial height configuration
        height[i,j] = h₀ + ϵ * randn()
    end
    th = zeros(size(height))
    th .= theta 
elseif device == "GPU"
    h = zeros(size(height))
    for i in 1:sys.Lx, j in 1:sys.Ly
        h[i,j] = h₀ + ϵ * randn()
    end
    # Lower contact angle inside the letters
	# Forward it to the GPU
	th = CUDA.adapt(CuArray, theta)
    height = CUDA.adapt(CuArray, h)
end
```

After the simulation box (or square to be more precise) is set we compute the first lattice Boltzmann equilibrium `Swalbe.equilibrium!(fout, height, velx, vely, vsq)`.
Knowing the initial equilibrium we can enter the lattice Boltzmann time loop.
Inside the loop we compute for every time step the forces that are present, here the film pressure (laplacian of the surface and wettability), slippage to regularize the contact line and the pressure gradient
```julia
Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, th, sys.n, sys.m, sys.hmin, sys.hcrit)
Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
# Forces are the pressure gradient and the slippage due to substrate liquid boundary conditions
Fx .= h∇px .+ slipx
Fy .= h∇py .+ slipy
```

Now that we know the forces we just have to update our distribution functions `fout` and `ftemp` (the hot sauce of the lattice Boltzmann method), in this case with a simple single relaxation time collision operator ([BGK](https://journals.aps.org/pr/abstract/10.1103/PhysRev.94.511)) and periodic boundary conditions. Last part of the lattice Boltzmann time step is the update of what is called macroscopic quantities (thickness & velocity), or simply the moment calculation (because these are the moments of the distribution mathematically speaking)
```julia
# Update the equilibrium
Swalbe.equilibrium!(feq, height, velx, vely, vsq)
# Collide and stream
Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
# New moments
Swalbe.moments!(height, velx, vely, fout)
```
and that's it.
Of course to generate data we make snapshots of the film using `Swalbe.snapshot!()` and return this collection of *thicknesses* at the end of the simulation.

What we get is something like this

![Hiern_logo_dewetting](https://user-images.githubusercontent.com/26249811/124448339-9cbc3880-dd82-11eb-9ccf-af44934b3f93.png)

All of the time steps that were generated during the simulation can be merged together and can be compressed into a movie, see below

![Dewetting_logo](https://gist.githubusercontent.com/Zitzeronion/807b9a7b2226e65643288df9a8cc1f46/raw/3a561e2a2b09eb42bf688f1d304f658b93fba8ed/logo_animation.gif)

This example will be further discussed in the Tutorials section.

## How to **perform research**

The numerical approach is quite robust for a lot of thin film simulations.
This means in the limit of small Reynolds (Re) and Mach number (Ma) simulations are usually stable.
As we employ the long wave approximation to derive the thin film equation we require that contact angle are small. 
Numerically the system is stable of up to 90$^{\circ}$ but to be within the theoretical assumptions we suggest to work with contact angles θ $\ll$ 1.

Things that we had a look at

- [Sliding droplets](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)
- [Surface tension gradients and droplet coalescence](https://arxiv.org/abs/2208.09307)
- Equilibrium shapes of droplets on patterned substrates, depending on surface tension.
- [Rayleigh-Taylor instability](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)
- Dewetting thin films
- [Dewetting of fluctuating thin films](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.034801)
- [Dewetting on switchable substrates](https://arxiv.org/abs/2112.09770)
- [Active thin films](http://dx.doi.org/10.1002/admi.202400835)
- [Ring-rivulets](https://arxiv.org/abs/2412.10230)
- [Multilayer films](https://arxiv.org/abs/2409.16659)

Things I have not yet looked into

- Non-Newtonian fluids
- Surfactants
- Particles
- Multicomponent/Multiphase

## Publications

- [Lattice Boltzmann simulations of the dynamics of thin liquid films](https://open.fau.de/items/bc35f5d4-b157-4aad-82be-9d23825cb37a/full)

Cumulative thesis for the grade Dr.-Ing. of Stefan Zitz.
About fifty pages of introduction to thin film dynamics, the shallow water system and numerical methods for thin films.
Some discussion on the [lattice Boltzmann method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) followed by a conclusion and four published articles (see below).   


- [Lattice Boltzmann method for thin-liquid-film hydrodynamics](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)

In this Phys. Rev. E publication we introduce the core concepts of **Swalbe**. We therefore discuss the [lattice Boltzmann method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) for shallow water flows and derive a matching condition with the [thin film equation](https://en.wikipedia.org/wiki/Thin-film_equation). Based on the derived model we present benchmark simulations and show that numerical results agree well with theoretical predictions, e.g., Cox-Voinov law.     

- [Lattice Boltzmann simulations of stochastic thin film dewetting](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.034801)

In the second Phys. Rev. E publication we discusses the addition of a stochastic force term as an effective model for the stochastic thin film equation, see [Grün et al.](https://link.springer.com/article/10.1007/s10955-006-9028-8).
This approach works surprisingly well and was validated against theory. Beyond the validation we showed how a film reacts to both a chemical pattern and thermal fluctuations. Again we were able to justify multiple findings with simple theoretical derivations.

- [Swalbe.jl: A lattice Boltzmann solver for thin film hydrodynamics](https://joss.theoj.org/papers/10.21105/joss.04312)

After the first two publications we cleaned up the repository and added this documentation as well as [test cases](https://github.com/Zitzeronion/Swalbe.jl/tree/master/test). 
During this process we became aware about the [Journal of Open Source Software](https://joss.theoj.org/) which we wholeheartly support. 
Scientific software is of value, however, it is often seen as byproduct and meeting the actual rigor requirement for a publication only discussing software is far from trivial. 
A JOSS publication was therefore a good place to show that **Swalbe.jl** is an open source tool and that it has a well definded scientific goal, namely thin liquid film solutions or the approximation of the thin film equation.   

- [Controlling the dewetting morphologies of thin liquid films by switchable substrates](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.8.L122001)

For this letter in Phys. Rev. Fluids we investigate a dewetting thin film on a "switchable" substrate. 
The idea is that based on external stimuli, such as light, local surface chemistry can be altered. 
We build a toy model that effects the local wettability both in time and space and found a non-dimensional number that explained the morphology changes that we observe and were able to motivate further results with theoretical considerations. 
Figure 1 and some videos have rendered from simulation data using [blender](https://www.blender.org/). 

- [Chemically Reactive Thin Films: Dynamics and Stability](http://dx.doi.org/10.1002/admi.202400835)

This model was mainly developed by [Tilman Richter](https://www.linkedin.com/in/tilman-richter-81a97b2a4/?originalSubdomain=de) and deals with active components that locally influences the surface tension. 
In contrast to a surfactant model this model allows for both a surface tension increase and decrease based on a parameter.
That said this generality made the model somewhat complex and the derivation of eigenmodes is far from trivial.

- [Instabilities of ring-rivulets: Impact of substrate wettability](https://arxiv.org/abs/2412.10230)

This paper has been accpeted in April 2025.
We study the dynamics of a ring-rivulet which is basically a liquid torus cut in half.
This initial condition has several interesting features and it is only a transient state, because it either retracts to single droplet or breaks up and forms multiple droplets.
Our work adds to the current understanding of this object by further considering wettability.
We show that a wettability gradient can turn a retracting ring-rivulet into one that breaks up and vise versa. 

### Preprints:

- [To coalesce or not to coalesce: Droplets and surface tension gradients](https://arxiv.org/abs/2208.09307)
- [Modeling of thin liquid films with arbitrary many layers](https://arxiv.org/abs/2409.16659)


## How to **support and contribute**

First of all leave a star if you like the idea of the project and/or the content of the package.
Second you can support the project by actively using it and raising [issues](https://github.com/Zitzeronion/Swalbe.jl/issues).
Help is always very welcome, if you want to contribute open a [**PR**](https://github.com/Zitzeronion/Swalbe.jl/pulls) or raise an [issue](https://github.com/Zitzeronion/Swalbe.jl/issues) with a feature request (and if possible with a way how to include it).
Feel free to DM me on [Bluesky](https://bsky.app/profile/zitzeronion.bsky.social) if you have questions, I try to answer them all timely.
