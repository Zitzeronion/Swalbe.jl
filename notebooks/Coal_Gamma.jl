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
	
    using Plots, Revise, DataFrames, FileIO, Swalbe, DataFramesMeta
end

# ╔═╡ bb534270-0e59-4c41-a825-fd6dc0fb4a7e
md"# Coalescence of sessile droplets

We numerically study the coalescence dynamics of two liquid droplets.
The droplets are placed on top of a rigid substrate with contact angle $\theta = 20^{\circ}$ and are barely in contact with each other.
Due to minimization of energy this system is not stable, instead the two droplets will coalesce into a single one.
This behavior can be explained following a self-similar ansatz for the thickness of the touching region.
With the single point $h_0(t_0)$, the so-called bridge height that defines the minimum in film thickness at the start of the experiment.
In the low contact angle regime one finds a solution given by

$h_0(t) \propto t^{2/3}.$

We're going to take this case as a starting point for the following experiments.
Instead of analyzing the outlined problem further, we introduce on top a surface tension gradient.
This gradient can for example be introduced with the addition of surfactants. 
The surfactant concentration is then e.g. limited to either one of the two droplets.

## Numerical setup

Simulations or numerical experiments are performed using the solver called *[Swalbe.jl](https://github.com/Zitzeronion/Swalbe.jl)*.
As of now these simulations are work in progress and surface tension gardients are not yet part of the master branch, therefore we use *[DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl)* to  call the local version of *Swalbe*.
Understand this as a note of caution, this is not a laboratory experiment!
For data analysis purposes, two more packages are added.
On the one hand *[Plots.jl](https://github.com/JuliaPlots/Plots.jl)* for basic plotting features and *[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)* for data storage and analysis.
First *DrWatson* is loaded to access a local version of *Swalbe*, and subsequently we load *Plots*,*DataFrames* and *FileIO* to write data to the disk."

# ╔═╡ 54427765-643f-44fe-84e1-c7c67b2cfe0d
md"## Surface tension

A lot of research has been dedicated towards the understanding of droplet coalescence.
However, as far the authors know the case with a variable surface tension along the coalescence axis has not been studied yet.
This case can arise for example assuming that one of the droplets has a non-vanishing surfactant concentration or one uses a well-defined heating source to locally heat the substrate.
Of course, we assume two things.
First being that the concentration gradient does not drive the dynamics.
Second being that the temperature distribution is fast as compared to the coalescence dynamics.
Both scenarios nevertheless introduce a surface tension gradient.

We are going to study the effect of this gradient using a simple one dimensional model.
The concrete choice for the spatially resolved surface tension is given by either one of three functions,

$\gamma^{lin.}(x) = \gamma_0\Bigg(1 - \epsilon \frac{x}{L}\Bigg),$

where γ₀ is a characteristic surface tension value, ϵ is a small deviation and $L$ is the size of the system. 
Besides this linear model, we use a model given by

$\gamma^{heavi}(x) = \begin{cases}\gamma_0\qquad\qquad\text{for x < L/2}\\ \gamma_0(1-\epsilon)\quad \text{else}\end{cases},$

with a Heaviside step function at the touching point.
A third option is given, using a tangent hyperbolicus to smoothen out the discontinuity of the Heaviside function.
Therefore, the function reads as

$γ^{smooth}(x) = γ₀[s(x) + (1-s(x))(1-ϵ)],$

with $s(x)$ being

$s(x) = \Bigg|1 - \left[\frac{1}{2} + \frac{1}{2}\tanh\left(\frac{x - a}{b}\right)\right]\Bigg|,$

where $a$ and $b$ define the intersection point and the smoothing width, respectively.

### System settings

For the following experiments, we use a one dimensional grid with $L=1024$ grid points.
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
	gamma_labels = Dict(1 => "default", 2 => "linear", 3 => "step", 4 => "tanh")
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

To do an experiment, we simply call a function that contains the *LBM* iterations.
This function will take as input arguments the spatially resolved surface tension γ(x).
Having a single function to run the experiments is rather convenient, as we simply can loop over if for further data.

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

# ╔═╡ 6332a336-fe15-4fb9-949b-c7d8ebc03176
md"To collect data, the only thing that is left to do is to run the function with the various surface tension fields we created.
For two of the four experiments, we know what should happen.
In case of a constant surface tension, the droplets will merge and the bridge height h₀ should grow as 

$h_0(t) \propto t^{2/3},$

according to references [[1](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/0DE6817B486607109F5A32018951B5C0/S002211209900662Xa.pdf/coalescence_of_liquid_drops.pdf),
[2](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/CB56EE5EDDAFF2E9503CFEB0D2110662/S0022112003004646a.pdf/div-class-title-inviscid-coalescence-of-drops-div.pdf),
[3](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.95.164503?casa_token=kdvf_DF7sK0AAAAA%3AN1Kw_LdxZMAfmGjHF5Cj39M7tugTp1frIl7HWV32tQw19igd43ha4Sk18H6fpJf2ZHruWHmQWK48), 
[4](https://arxiv.org/pdf/1406.5842.pdf)]. 

The other extreme is the step function.
One can create this simply by using two different but miscible liquids for the different droplets. 
Naturally, the different liquids admit different surface tensions, and therefor there is a step like profile in the surface tension.
This case has been studied by Karpitschka et al. [[5](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/sharp-transition-between-coalescence-and-noncoalescence-of-sessile-drops/121D63B9ABFD9040C327ED8AC10FBE01),
[6](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.066103)] and relates to interesting result

$h_0(t) \sim h_0(t_0),$

thus the bridge is not growing but keeps the droplet separated.
While the capillarity tries to merge the droplets to reduce the overall curvature due to the neck, a Marangoni like flow counteracts the capillarity.
This counteracting forces stabilizes the two droplet state.

We assume that the other two surface tension gradients, both the linear and the smoothed step function, will lead to coalescence. 
However, especially in the case of smoothing function using a `tanh` should be dependent on the smoothing width.
But first we have to perform the experiments.

## Numerical experiments
"

# ╔═╡ a0b3c869-3a7e-4b10-a2d8-7021b8c4c54d
function plot_data(data; k=1, t1=100, t2=500, t3=1000, leg=true)
	tarray = 100:100:1000000000
	if isa(data, Array) 
		plot(data[k, t1, :], 
			label="γ=$(gamma_labels[k]) t=$(tarray[t1])", 
			line = (:auto, 4), 
			xlabel="x", 
			ylabel="h",
			legendfontsize=14,
			guidefont = (18, :black), 
			tickfont = (12, :black),
			legend= leg,
			xlims=(472, 552)
		)
		for i in [t2, t3]
			plot!(data[k, i, :], label="γ=$(gamma_labels[k]) t=$(tarray[i])", line = (:auto, 4))
		end
		ylims!(0,15)
	elseif isa(data, Swalbe.SysConst_1D)
		save_file = "..\\data\\Drop_coalescence\\gamma_$(gamma_labels[k])_tmax_$(data.param.Tmax).jld2"
		df = load(save_file) |> DataFrame
		plot(df[!, Symbol("h_$(tarray[t1])")], 
			label="γ=$(gamma_labels[k]) t=$(tarray[t1])", 
			line = (:auto, 4), 
			xlabel="x", 
			ylabel="h",
			legendfontsize=14,
			guidefont = (18, :black), 
			tickfont = (12, :black),
			legend= leg,
			xlims=(472, 552)
		)
		for i in [t2, t3]
			plot!(df[!, Symbol("h_$(tarray[i])")], label="γ=$(gamma_labels[k]) t=$(tarray[i])", line = (:auto, 4))
		end
		ylims!(0,15)
	end
end

# ╔═╡ 2f6e1154-eaff-4c20-9fac-290474f45f0b
md"### First sweep

We let the system evolve up to a 10⁵ time steps.
What we hope to observe is that the bridge does not grow for the step like surface tension gradient.

"

# ╔═╡ 2edc58c6-4ee0-4c5e-8013-311e81820c4c
begin
	sys = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=100000, δ=10.0))
	
	data = zeros(4, 1000, 1024)
	
	for i in 1:4
	 	data[i, :, :] = run_(sys, γ[i, :], r₁=rad, r₂=rad)

		println("Done with iteration $i")
	end
end

# ╔═╡ bab2a9ab-1ee1-4146-95e7-03d87a8f9c35
md"Hard to see, but the lower left plot is the one generated with step like surface tension.
Indeed, the bridge does not grow, but interestingly the position of the neck is moving.
It is moving in the direction of higher surface tension, in fact.

The other three panels show a somewhat similar evolution.
Somewhat because the evolution on the right side is not fully symmetric around the center.
But we will come back to this."

# ╔═╡ 35541b43-a834-4f01-ad0d-fa11be9af74b
begin
	p11=plot_data(data, k=1)
	p12=plot_data(data, k=2)
	p13=plot_data(data, k=3)
	p14=plot_data(data, k=4)
	plot(p11, p12, p13, p14, legend = false)
end

# ╔═╡ b164e7ec-eaa8-4f2a-b35a-6ac01fd12875
md"### Second sweep

Now that we are, at least to some degree, happy with the short time evolution of the bridges and thus the (non)coalescence process, we move on to collect more data, especially towards later stages of the coalescence process.

This is what we do in the second sweep. 
Instead of stopping after 10⁵ time steps, we take up 5x10⁶ time steps into consideration.
Okay but what does a single time step mean?
The short answer to that question is, it means nothing.
To make sense of this time steps we use characteristic time scales and match them with our parameters. 
We proceed in a similar matter for the lenght scale. 
thus we have $\tau$ and $l$ with defining equations,

$\tau = \frac{r\mu}{\gamma},$

where $\tau$ is the time scale of viscos relaxation and thus $t/\tau$ is nondimensional.
For the length scale we get quite naturally

$l = r,$

where $r$ is the base radius of the droplet.
The numerical values therefore are

| Parameter   |      Value      | 
|:----------:|:-------------:|
| μ |  1/6 |
| γ₀ | 10⁻⁴ |
| l | 171 |
| τ | 285000 | 

The value of $\tau$ can now be used to make sense of our time steps.
In our second sweep we simulate for 5×10⁶ time steps and in terms of nondimensional time up to 17.5[t/τ].
A single time step therefor creates the increment of 3.5×10⁻⁶ in terms of t/τ.
Which is quite cool because for scaling laws we can cover a huge range in time, from 10⁻⁵ up to 10.
"

# ╔═╡ 9e714346-0e8d-4de5-85b4-4ede3e275834
"""
	tau_vr(;r₀=171, μ=1/6, γ=γ₀)

Computation of the viscous relaxation time.
"""
function tau_vr(;r₀=171, μ=1/6, γ=γ₀)
	return r₀ * μ / γ
end

# ╔═╡ a2389b72-c460-4026-8676-b3b0fb4acdc2
"""
	tau_ic(;ρ=1, r₀=171, γ=γ₀)

Computation of the inertio-capillary time.
"""
function tau_ic(;ρ=1, r₀=171, γ=γ₀)
	return sqrt(ρ*r₀^3/γ)
end

# ╔═╡ 9e567d24-9636-467a-b491-ed9bdc584854
begin
	sys2 = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=5000000, δ=1.0))
	# Memory
	data2 = zeros(4, 50000, 1024)
	df2 = DataFrame()
	# Time
	time = 100:100:size(data2)[2]*100
	# Loop over γ
	for i in 1:4
		# Check if there is already a file created
		save_file = "..\\data\\Drop_coalescence\\gamma_$(gamma_labels[i])_tmax_$(sys2.param.Tmax).jld2"
		# If so just read it from disc
		if isfile(save_file)
			tmp = load(save_file) |> DataFrame
			tmp.nablaG .= gamma_labels[i]
			df2 = vcat(df2, tmp)
		# If not, compute the evolution of the droplet coalescence
		else
	 		data2[i, :, :] = run_(sys2, γ[i, :], r₁=rad, r₂=rad)
			# Think about storing them on the disc
			df_fluid = Dict()
			# Loop through the time once more
			for t in 1:size(data2)[2]
		        df_fluid["h_$(time[t])"] = data2[i,t,:]
		    end
			# Add a surface tension field label
			df_fluid["nablaG"] = fill(gamma_labels[i], 1024)
			# Save it
			println("Writing file $(gamma_labels[i]) to disk")
			save("..\\data\\Drop_coalescence\\gamma_$(gamma_labels[i])_tmax_$(sys2.param.Tmax).jld2", df_fluid)
		end
		# Print that you are done
		println("Done with iteration $(gamma_labels[i])")
	end
end

# ╔═╡ 19772481-02f9-4091-bfc9-e2853e64d4d6
md"#### Results - Inspection

Two things come to mind imideatly, first being the now pronounced asymmetry around the center for the experiments on right.
In the lower right this asymmetry manifest itself already at the time step of the orange curve.

The second eyecatcher is the evolution of the lower left.
Similar to the short experiments the neck is not growing, but moves towards the region of higher surface tension.
Therefor, the trend of noncoalescing with the step function surface tension is in fact a stable one.
The two droplets keep apart for the entirety of the experiments.
"

# ╔═╡ 0d7c9262-ff25-46bd-9833-bddfd7f958f9
begin
	p1=plot_data(sys2, k=1, t2=5000, t3=45000)
	p2=plot_data(sys2, k=2, t2=5000, t3=45000)
	p3=plot_data(sys2, k=3, t2=5000, t3=45000)
	p4=plot_data(sys2, k=4, t2=5000, t3=45000)
	plot(p1, p2, p3, p4, legend = false)
end

# ╔═╡ 50cb438d-501e-411e-844d-e75569ca3c85
md"#### Results - Bridge height

After a first visual inspection of the neck region we move towards a more quantitative analysis.
We would like to know how the height of the neck is evolving with time.
This is quite a simple measurement and will tell us already a lot about our system.

"

# ╔═╡ e5b058aa-16ca-4864-84d0-9879b353aa18
@subset(df2, :nablaG .== "tanh")

# ╔═╡ b8eadb42-643c-41e3-ae94-7fe0cefb3d6b
function bridge_height(df; time=100:100:5000000)
	df_ = DataFrame()
	time_list = []
	height_list = []
	grad_list = []
	for i in gamma_labels
		tmp = @subset(df2, :nablaG .== i) 
		for t in time
			push!(time_list, t)
			push!(grad_list, i)
			push!(height_list, minimum(tmp[!, Symbol("h_$(t)")][]))
end

# ╔═╡ d591ec33-0ae9-4d64-9422-9a1f283de7b5
md"### Data

Instead of rerunning the experiments, we write them to disc and relaod them.
While it does not save a lot of time, given the experiments run not even for ten minutes, we have the ability to share the raw data if needed."

# ╔═╡ f5696349-24f8-4503-b037-ad511579918d
# Save data, so it needs to be loaded not computed

# ╔═╡ 3091f849-7ce2-400d-9938-4b89215a0bdc
md"## References

1. Eggers, J., Lister, J.R. and Stone, H.A., 1999. Coalescence of liquid drops. *Journal of Fluid Mechanics*, *401*, pp.293-310.
2. Duchemin, L., Eggers, J. and Josserand, C., 2003. Inviscid coalescence of drops. *Journal of Fluid Mechanics*, *487*, pp.167-178.
3. Aarts, D.G., Lekkerkerker, H.N., Guo, H., Wegdam, G.H. and Bonn, D., 2005. Hydrodynamics of droplet coalescence. *Physical review letters*, *95*(16), p.164503.
4. Sprittles, J.E. and Shikhmurzaev, Y.D., 2014. A parametric study of the coalescence of liquid drops in a viscous gas. *Journal of Fluid Mechanics*, *753*, pp.279-306.
5. Karpitschka, S., & Riegler, H., 2014. Sharp transition between coalescence and non-coalescence of sessile drops. *Journal of Fluid Mechanics*, 743.
6. Karpitschka, S. and Riegler, H., 2012. Noncoalescence of sessile drops from different but miscible liquids: hydrodynamic analysis of the twin drop contour as a self-stabilizing traveling wave. *Physical review letters*, *109*(6), p.066103..
7. Eddi, A., Winkels, K.G. and Snoeijer, J.H., 2013. Influence of droplet geometry on the coalescence of low viscosity drops. *Physical review letters*, *111*(14), p.144502.
8. Zitz, S., Scagliarini, A., Maddu, S., Darhuber, A.A. and Harting, J., 2019. Lattice Boltzmann method for thin-liquid-film hydrodynamics. *Physical Review E*, *100*(3), p.033313.
9. Oron, A., Davis, S.H. and Bankoff, S.G., 1997. Long-scale evolution of thin liquid films. Reviews of modern physics, 69(3), p.931."

# ╔═╡ Cell order:
# ╟─bb534270-0e59-4c41-a825-fd6dc0fb4a7e
# ╠═52725098-857a-4301-b86b-d9cd819de541
# ╟─54427765-643f-44fe-84e1-c7c67b2cfe0d
# ╠═eda3dc93-7626-42eb-82a6-b8615bd0f477
# ╠═bcb187be-9a59-46ce-aebe-74e7003077d8
# ╠═677ff3cc-4037-4b19-a521-dbca74a635a7
# ╟─c4236e13-fca3-4350-adf7-b98c7bde8a0a
# ╟─3594c0d9-0010-4086-9e7e-163bdf1b0195
# ╟─708f54fc-0bd4-4577-85ec-4faf38029c2f
# ╟─8010c641-a385-4f3f-a88d-817332e45091
# ╠═09a80dac-0cd5-42f3-9676-e412a58f58db
# ╟─ac41b37e-841f-47c6-b5ff-3b10fc2c86ae
# ╠═547e2ffb-b0a9-4bf2-a80a-5a6b5aed7e5a
# ╟─6332a336-fe15-4fb9-949b-c7d8ebc03176
# ╟─a0b3c869-3a7e-4b10-a2d8-7021b8c4c54d
# ╟─2f6e1154-eaff-4c20-9fac-290474f45f0b
# ╠═2edc58c6-4ee0-4c5e-8013-311e81820c4c
# ╟─bab2a9ab-1ee1-4146-95e7-03d87a8f9c35
# ╠═35541b43-a834-4f01-ad0d-fa11be9af74b
# ╟─b164e7ec-eaa8-4f2a-b35a-6ac01fd12875
# ╟─9e714346-0e8d-4de5-85b4-4ede3e275834
# ╟─a2389b72-c460-4026-8676-b3b0fb4acdc2
# ╠═9e567d24-9636-467a-b491-ed9bdc584854
# ╟─19772481-02f9-4091-bfc9-e2853e64d4d6
# ╠═0d7c9262-ff25-46bd-9833-bddfd7f958f9
# ╠═50cb438d-501e-411e-844d-e75569ca3c85
# ╠═e5b058aa-16ca-4864-84d0-9879b353aa18
# ╠═b8eadb42-643c-41e3-ae94-7fe0cefb3d6b
# ╟─d591ec33-0ae9-4d64-9422-9a1f283de7b5
# ╠═f5696349-24f8-4503-b037-ad511579918d
# ╟─3091f849-7ce2-400d-9938-4b89215a0bdc
