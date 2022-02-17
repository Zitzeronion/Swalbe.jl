### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 52725098-857a-4301-b86b-d9cd819de541
using DrWatson

# ╔═╡ 005b7402-df13-410c-88e1-69fec65b5ff0
# Because I am working on a local branch
@quickactivate :Swalbe

# ╔═╡ f75e02be-8f28-11ec-16ae-8bc62cd08643
using Plots, DataFrames, FileIO

# ╔═╡ bb534270-0e59-4c41-a825-fd6dc0fb4a7e
md"# Coalescence of sessile droplets

We numerically study the coalescence dynamics of two liquid droplets.
The droplets are placed on top of a rigid substrate with contact angle $\theta = 20^{\circ}$ and are barly in contact with each other.
Due to minimization of energy this system is not stabel, instead the two droplets will coalesce into a single one.
This behavior can be explained following a self-similar ansatz for the thickness of the touching region.
With the single point $h_0(t_0)$, the socalled bridge height that defines the minimum in film thickness at the start of the experiment.
In the low contact angle regime one finds a solution givenn by

$h_0(t) \propto t^{2/3}.$

We gona take this case as starting point for the following experiments.
Instead of analyzing the outlined problem further we introduce on top a surface tension gradient.
This gradient can for example be introduce with the addition of surfactants. 
The surfactant concentration is then e.g. limited to either one of the two droplets.

## Numerical setup

Simulations or numerical experiments are performed using the solver called *[Swalbe.jl](https://github.com/Zitzeronion/Swalbe.jl)*.
As of now this simulations are work in progress and surface tension gardients are not yet part of the master branch, therefore we use *[DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl)* to  call the local version of *Swalbe*.
Understand this as a note of caution, this is not a laboratory experiment!
For data analysis purposes two more packages are added.
On the one hand *[Plots.jl](https://github.com/JuliaPlots/Plots.jl)* for basic plotting features and *[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)* for data storage and analysis.
First *DrWatson* is loaded to access a local version of *Swalbe* and subsequently we load *Plots*,*DataFrames* and *FileIO* to write data to the disk."

# ╔═╡ 54427765-643f-44fe-84e1-c7c67b2cfe0d
md"## Surface tension

A lot of research has been dedicated towards the understanding of droplet coalescence.
However as far the authors know the case with a variable surface tension along the coalescence axis has not been studied yet.
This case can arise for example assuming that one of the droplets has a non vanishing surfactant concentration or one uses a well defined heating source to locally heat the substrate.
Of course we assume two things.
First being that the concentration gradient does not drive the dynamics.
Second being that the temperature distribution is fast as compared to the coalescence dynamics.
Both scenarios nevertheless introduce a surface tension gardient.

We are going to study the effect of this gradient using a simple one dimensional model.
The concrete choice for the spatially resolved surface tension is given by either one of three functions,

$\gamma^{lin.}(x) = \gamma_0\Bigg(1 - \epsilon \frac{x}{L}\Bigg),$

where γ₀ is a characteristic surface tension value, ϵ is a small deviation and $L$ is the size of the system. 
Besides this linear model we us a model given by

$\gamma^{heavi}(x) = \begin{cases}\gamma_0\qquad\qquad\text{for x < L/2}\\ \gamma_0(1-\epsilon)\quad \text{else}\end{cases},$

with a Heaviside step function at the touching point.
A third option is given using a tagent hyperbolicus to smoothen out the discontinuity of the Heaviside function.
Therefore the function reads as

$γ^{smooth}(x) = γ₀[s(x) + (1-s(x))(1-ϵ)],$

with $s(x)$ being

$s(x) = \Bigg|1 - \left[\frac{1}{2} + \frac{1}{2}\tanh\left(\frac{x - a}{b}\right)\right]\Bigg|,$

where $a$ and $b$ define the intersection point and the smoothing width respectively.

### System settings

For the following experiments we use a one dimensional grid with $L=1024$ grid points.
We store the three $\gamma$ functions in a single array and will loop through them. 
"

# ╔═╡ eda3dc93-7626-42eb-82a6-b8615bd0f477
begin
	L = 1024
	x = zeros(L)
	γ = zeros(3,L)
	γ₀ = 0.0001
	smooth = abs.(1 .- (0.5 .+ 0.5 .* tanh.((collect(1:L) .- L÷2) ./ (L÷10))))
end

# ╔═╡ bcb187be-9a59-46ce-aebe-74e7003077d8
function gamma_curves!(x; x0=γ₀, ϵ=0.1)
	x[1,:] .= x0 .* (1 .- ϵ .* collect(1:L) ./ L)
	x[2,1:L÷2] .= x0
	x[2,L÷2+1:L] .= x0 - x0 * ϵ
	# x[3,:] .= x0 .* exp.(-collect(1:L)/L)
	x[3,:] .= x0 .* smooth .+ (1 .- smooth) .* x0 .*(1 - ϵ) 
end

# ╔═╡ 677ff3cc-4037-4b19-a521-dbca74a635a7
gamma_curves!(γ)

# ╔═╡ c4236e13-fca3-4350-adf7-b98c7bde8a0a
@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 10)
    n = length(y)
    sx, sy = x[1:step:n], y[1:step:n]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := [Inf]
        y := [Inf]
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end

# ╔═╡ 3594c0d9-0010-4086-9e7e-163bdf1b0195
md"Below is a plot of the three choices of $\gamma(x)$ normalized with the characteristic value $\gamma_0$. "

# ╔═╡ 708f54fc-0bd4-4577-85ec-4faf38029c2f
begin
	smap = 50
	plot(collect(1:L), γ[1,:] ./ γ₀, 
		 w=3, 
		 st = :samplemarkers,
		 step = smap, 						
		 marker = (8, :auto, 0.6),			
		 label="lin", 
		 xlabel="x/Δx", 
		 ylabel="γ/γ₀", 
	     legendfontsize = 14,			# legend font size
         tickfont = (14),	# tick font and size
         guidefont = (15)	# label font and size
		 )
	plot!(collect(1:L), γ[2,:] ./ γ₀, 
		  w=3, 
		  label="step",
		  st = :samplemarkers,
		  step = smap, 						
		  marker = (8, :auto, 0.6),)
	plot!(collect(1:L), γ[3,:] ./ γ₀, 
		  w=3, 
		  label="tanh",
	      st = :samplemarkers,
		  step = smap, 						
		  marker = (8, :auto, 0.6))
end

# ╔═╡ 8010c641-a385-4f3f-a88d-817332e45091
md"We then apply this surface tension function $\gamma(x)$ to the initial state displayed below.
The thickness of the two droplets is given by $h(x)$, with a touching point in the center of the domain at $L/2$."

# ╔═╡ 09a80dac-0cd5-42f3-9676-e412a58f58db
	begin
	# The initial configuration of the numerical experiment
	rad = 500
	h = Swalbe.two_droplets(Swalbe.SysConst_1D(L=1024), r₁=rad, r₂=rad)
	p0 = plot(collect(1:1024), h, 
		      w=3, 
		      aspect_ratio=7, 
		      label="Initial conf.", 
		      xlabel="x [Δx]", 
		      ylabel="h(x)",
		      legendfontsize = 14,			# legend font size
              tickfont = (14),	# tick font and size
              guidefont = (15)	# label font and size
	          )
	ylims!(0,40)
end

# ╔═╡ ac41b37e-841f-47c6-b5ff-3b10fc2c86ae
md"### Run function

To do an experiment we simply call a function that encapsules the *LBM* iterations.
This function will take as input arguments the surface tension as well as the pair $(n,m)$ for the disjoining pressure $\Pi(h)$.
Due to the fact that these are input variables can simply loop over them and collect all the data we are interested in.

Below is the definition of the function `run_()`"

# ╔═╡ 547e2ffb-b0a9-4bf2-a80a-5a6b5aed7e5a
function run_(
    sys::Swalbe.SysConst_1D,
    gamma::Vector;
    r₁=115,
    r₂=115, 
    θ₀=1/9,  
    verbos=true, 
    dump = 100, 
    fluid=zeros(sys.Tmax÷dump, sys.L)
)
    println("Simulating droplet coalecense with surface tension gardient")
    state = Swalbe.Sys(sys, kind="gamma")
    drop_cent = (sys.L/3, 2*sys.L/3)
    state.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state)
    state.γ .= gamma
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            if verbos
                println("Time step $t bridge height is $(round(state.height[Int(sys.L/2)], digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(fluid, state.height, t, dumping = dump)
    end

    return fluid
    
end

# ╔═╡ 2edc58c6-4ee0-4c5e-8013-311e81820c4c
begin
	sys = Swalbe.SysConst_1D(L=1024, Tmax=4000000, δ=50.0)
	l = run_(sys, fill(1e-4, 1024), r₁=rad, r₂=rad)
end

# ╔═╡ c2753149-8036-49d6-86e0-ea3d8fccf7cf
begin
	# Sphere radius
	# sphere_rad = 500
	# Different surface tension values
	# γs = [0.00008, 0.0003, 0.0006, 0.0008]
	# Array that contains the simulation data
	# data_merge = zeros(20000, 1024, length(γs))
	# Loop different surface tension values
	#for γ in enumerate(γs)
		# System parameter, δ=50 can still considered small to medium slippage
		# sys = Swalbe.SysConst_1D(L=1024, Tmax=4000000, δ=50.0)
		# The experiment
		# data_merge[:,:,γ[1]] = run_drop_coal(sys, r₁=sphere_rad, r₂=sphere_rad)
	# end
end

# ╔═╡ Cell order:
# ╟─bb534270-0e59-4c41-a825-fd6dc0fb4a7e
# ╠═52725098-857a-4301-b86b-d9cd819de541
# ╠═005b7402-df13-410c-88e1-69fec65b5ff0
# ╠═f75e02be-8f28-11ec-16ae-8bc62cd08643
# ╟─54427765-643f-44fe-84e1-c7c67b2cfe0d
# ╠═eda3dc93-7626-42eb-82a6-b8615bd0f477
# ╠═bcb187be-9a59-46ce-aebe-74e7003077d8
# ╠═677ff3cc-4037-4b19-a521-dbca74a635a7
# ╟─c4236e13-fca3-4350-adf7-b98c7bde8a0a
# ╟─3594c0d9-0010-4086-9e7e-163bdf1b0195
# ╟─708f54fc-0bd4-4577-85ec-4faf38029c2f
# ╟─8010c641-a385-4f3f-a88d-817332e45091
# ╠═09a80dac-0cd5-42f3-9676-e412a58f58db
# ╠═ac41b37e-841f-47c6-b5ff-3b10fc2c86ae
# ╠═547e2ffb-b0a9-4bf2-a80a-5a6b5aed7e5a
# ╠═2edc58c6-4ee0-4c5e-8013-311e81820c4c
# ╠═c2753149-8036-49d6-86e0-ea3d8fccf7cf
