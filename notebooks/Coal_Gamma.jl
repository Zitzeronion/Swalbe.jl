### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 52725098-857a-4301-b86b-d9cd819de541
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	
    using Plots, Revise, DataFrames, FileIO, Swalbe
end

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
	x = collect(1:L)
	γ = zeros(4,L)
	ε = 0.2
	γ₀ = 0.0001
	γ_bar = (γ₀ + (γ₀ - ε))/2
	Δγ = ε
	sl = L÷10
end

# ╔═╡ bcb187be-9a59-46ce-aebe-74e7003077d8
function gamma_curves!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	out[1,:] .= x0
	out[2,:] .= x0 .* (1 .- ϵ .* l ./ L)
	out[3,1:L÷2] .= x0
	out[3,L÷2+1:L] .= x0 - x0 * ϵ
	# x[3,:] .= x0 .* exp.(-collect(1:L)/L)
	out[4,:] .= x0 .* smooth(x, L, sl) .+ (1 .- smooth(x, L, sl)) .* x0 .*(1 - ϵ) 
	return nothing
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
	plot(collect(1:L), γ[2,:] ./ γ₀, 
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
	plot!(collect(1:L), γ[3,:] ./ γ₀, 
		  w=3, 
		  label="step",
		  st = :samplemarkers,
		  step = smap, 						
		  marker = (8, :auto, 0.6),)
	plot!(collect(1:L), γ[4,:] ./ γ₀, 
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
	h = Swalbe.two_droplets(Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs()), r₁=rad, r₂=rad)
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
    fluid=zeros(sys.param.Tmax÷dump, sys.L)
)
    println("Simulating droplet coalecense with surface tension gardient")
    state = Swalbe.Sys(sys, kind="gamma")
    drop_cent = (sys.L/3, 2*sys.L/3)
    state.basestate.height .= Swalbe.two_droplets(sys, r₁=r₁, r₂=r₂, θ₁=θ₀, θ₂=θ₀, center=drop_cent)
    Swalbe.equilibrium!(state, sys)
    state.γ .= gamma
	Swalbe.∇γ!(state)
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.param.Tmax
        if t % sys.param.tdump == 0
            mass = 0.0
            mass = sum(state.basestate.height)
            if verbos
                println("Time step $t bridge height is $(round(minimum(state.basestate.height[sys.L÷2-20:sys.L÷2+20]), digits=3))")
            end
        end
        Swalbe.filmpressure!(state, sys, γ=gamma)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.basestate.F .= -state.basestate.h∇p .- state.basestate.slip .- state.∇γ
        Swalbe.equilibrium!(state, sys)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        
        Swalbe.snapshot!(fluid, state.basestate.height, t, dumping = dump)
    end
	# println("Max γ: $(maximum(state.γ)),\nMin γ: $(minimum(state.γ))")
    return fluid
    
end

# ╔═╡ f0f554f6-d9ab-4eba-a186-a07f481904cb
sys = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=100000, δ=10.0))

# ╔═╡ 2edc58c6-4ee0-4c5e-8013-311e81820c4c
begin
	data = zeros(4, 1000, 1024)
	for i in 1:4
	 	data[i, :, :] = run_(sys, γ[i, :], r₁=rad, r₂=rad)
		println("Done with iteration $i")
	end
end

# ╔═╡ a0b3c869-3a7e-4b10-a2d8-7021b8c4c54d
function plot_data(data, k, t1, t2, t3)
	plot(data[k, t1, :], label="γ=lin. t=$(t1*100)", line = (:auto, 4), xlabel="x", ylabel="h", xlims=(472, 552))
	for i in [t2, t3]
		plot!(data[k, i, :], label="γ=lin. t=$(i*100)", line = (:auto, 4))
	end
	ylims!(0,15)
end

# ╔═╡ c5118d35-2015-49ee-889a-2e3040e906eb
begin
	t1 = 100
	t2 = 500
	t3 = 1000
	
	plot_data(data, 1, t1, t2, t3)

end

# ╔═╡ 6922371e-4418-46ae-9f39-7690f78e8b45
plot_data(data, 3, t1, t2, t3)

# ╔═╡ f992a3c6-8eca-47e6-bedb-1a2486b1a04e
plot_data(data, 4, t1, t2, t3)

# ╔═╡ Cell order:
# ╟─bb534270-0e59-4c41-a825-fd6dc0fb4a7e
# ╠═52725098-857a-4301-b86b-d9cd819de541
# ╟─54427765-643f-44fe-84e1-c7c67b2cfe0d
# ╠═eda3dc93-7626-42eb-82a6-b8615bd0f477
# ╠═bcb187be-9a59-46ce-aebe-74e7003077d8
# ╠═677ff3cc-4037-4b19-a521-dbca74a635a7
# ╟─c4236e13-fca3-4350-adf7-b98c7bde8a0a
# ╟─3594c0d9-0010-4086-9e7e-163bdf1b0195
# ╠═708f54fc-0bd4-4577-85ec-4faf38029c2f
# ╟─8010c641-a385-4f3f-a88d-817332e45091
# ╠═09a80dac-0cd5-42f3-9676-e412a58f58db
# ╟─ac41b37e-841f-47c6-b5ff-3b10fc2c86ae
# ╠═547e2ffb-b0a9-4bf2-a80a-5a6b5aed7e5a
# ╠═f0f554f6-d9ab-4eba-a186-a07f481904cb
# ╠═2edc58c6-4ee0-4c5e-8013-311e81820c4c
# ╠═a0b3c869-3a7e-4b10-a2d8-7021b8c4c54d
# ╠═c5118d35-2015-49ee-889a-2e3040e906eb
# ╠═6922371e-4418-46ae-9f39-7690f78e8b45
# ╠═f992a3c6-8eca-47e6-bedb-1a2486b1a04e
