### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ 52725098-857a-4301-b86b-d9cd819de541
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	
    using Swalbe
end

# ╔═╡ 691876ad-2ee6-4b87-974f-66a3650c4b2f
using Plots, Revise, DataFrames, FileIO, DataFramesMeta, StatsBase, CSV

# ╔═╡ eda3dc93-7626-42eb-82a6-b8615bd0f477
include("helpers.jl")

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
         tickfontsize = 14,	# tick font and size
         guidefontsize = 15	# label font and size
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

# ╔═╡ 28025793-b001-4597-aa1c-f2dd06c8a34e
begin
	# Some parameters and the initial condition
	rad = 500
	h = Swalbe.two_droplets(Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs()), r₁=rad, r₂=rad)
end

# ╔═╡ 76aff625-dcfb-4875-8366-d4c6aac51e54
begin
	# The initial configuration of the numerical experiment
	init_p = plot(collect(1:1024), h, 
		      w=3, 
		      aspect_ratio=5, 
		      label="", 
			  title="Initial condition",
			  ribbon=(h,0),
		      xlabel="x [Δx]", 
		      ylabel="h(x) [l.b.u.]",
		      legendfontsize = 14,			# legend font size
              tickfontsize = 14,	# tick font and size
              guidefontsize = 15,	# label font and size
			  grid = false,
	          )
	ylims!(0,40)
	xlims!(0,1024)
end

# ╔═╡ 77566ed2-14b6-48c3-831f-80efce8c2e3e
# savefig(init_p, "..\\..\\figures\\initial_cond.svg")

# ╔═╡ ac41b37e-841f-47c6-b5ff-3b10fc2c86ae
md"### Run function

To do an experiment, we simply call a function that contains the *LBM* iterations.
This function will take as input arguments the spatially resolved surface tension γ(x).
Having a single function to run the experiments is rather convenient. 
This way we simply can loop over for further data collection.

The definition of the function can be found in Swalbe.jl (*surface\_tension\_gradient* branch) and is called `run_gamma()`

To collect data, the only thing that is left to do is to run the function with the various surface tension fields we created.
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

*Side note*, the experiments presented by Karpitschka et al. were verified using a numerical solver for the thin film problem.
They numerically solved 

$\partial_t h = -\partial_x\bigg\{\frac{1}{\mu}\bigg[\frac{h^3}{3}\partial_x(\gamma\partial_x^2h)+\frac{h^2}{2}\partial_x\gamma\bigg]\bigg\},$

with $\mu$ being the viscosity.
This model is more than appropriate for macroscopic drops.
However if the overall size of the problem gets smaller another term could contribute, namely the disjoining pressure,

$\Pi(h) = \kappa(\theta)\bigg[\bigg(\frac{h_{\ast}}{h}\bigg)^n - \bigg(\frac{h_{\ast}}{h}\bigg)^m\bigg],$

where $\kappa(\theta)$ is a function of the contact angle $\theta$, $h_{\ast}$ a parameter for which $\Pi(h_{\ast}) = 0$ and $n, m$ are positive integers obeying $n > m > 0$.
Combining this term and the second derivative of the thickness into single expression 

$p_{\text{film}} = \gamma \partial_x^2 h + \Pi(h),$

the thin film equation becomes

$\partial_t h = -\partial_x\bigg\{\frac{1}{\mu}\bigg[\frac{h^3}{3}\partial_x p +\frac{h^2}{2}\partial_x\gamma\bigg]\bigg\},$

while this may not seem like a huge game changer, the disjoining pressure is there to actually account for things like the three phase contact line.
We shall find out now if we are able to come to the same conclusion as Karpitschka with the addition of the disjoining pressure.

## Numerical experiments
"

# ╔═╡ a0b3c869-3a7e-4b10-a2d8-7021b8c4c54d
function plot_data(data; k=1, t1=100, t2=500, t3=1000, leg=true, labels=gamma_labels)
	tarray = 100:100:1000000000
	if isa(data, Array) 
		plot(data[k, t1, :], 
			label="γ=$(labels[k]) t=$(tarray[t1])", 
			line = (:auto, 4), 
			xlabel="x", 
			ylabel="h",
			legendfontsize=14,
			guidefontsize = 18, 
			tickfontsize = 12,
			legend= leg,
			xlims=(472, 552)
		)
		for i in [t2, t3]
			plot!(data[k, i, :], label="γ=$(labels[k]) t=$(tarray[i])", line = (:auto, 4))
		end
		ylims!(0,15)
	elseif isa(data, Swalbe.SysConst_1D)
		save_file = "..\\..\\data\\Drop_coalescence\\gamma_$(labels[k])_tmax_$(data.param.Tmax).jld2"
		df = load(save_file) |> DataFrame
		plot(df[!, Symbol("h_$(tarray[t1])")], 
			label="γ=$(labels[k]) t=$(tarray[t1])", 
			line = (:auto, 4), 
			xlabel="x", 
			ylabel="h",
			legendfontsize=14,
			guidefontsize = 18, 
			tickfontsize = 12,
			legend= leg,
			xlims=(472, 552)
		)
		for i in [t2, t3]
			plot!(df[!, Symbol("h_$(tarray[i])")], label="γ=$(labels[k]) t=$(tarray[i])", line = (:auto, 4))
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
	 	data[i, :, :] = Swalbe.run_gamma(sys, γ[i, :], r₁=rad, r₂=rad)

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

$\tau_v = \frac{l\mu}{\gamma},\qquad \tau_i = \sqrt{\frac{\rho l^3}{\gamma}},$

where $\tau_v$ is the time scale of viscos relaxation and $\tau_i$ is the inertio-capillary time.
The product $t/\tau$ is therefore nondimensional.
For the length scale we get quite naturally

$l = r,$

where $r$ is the base radius of the droplet.
The numerical values therefore are

| Parameter   |      Value      | 
|:----------:|:-------------:|
| μ |  1/6 |
| γ₀ | 10⁻⁴ |
| l | 171 |
| τᵥ | 285000 |
| τᵢ | 223611.5 |

The value of $\tau$ (because $\tau_i \approx \tau_v$) can now be used to make sense of our time steps.
In our second sweep we simulate for 5×10⁶ time steps and in terms of nondimensional time up to 17.5[t/τ].
A single time step therefor creates the increment of 3.5×10⁻⁶ in terms of t/τ.
Which is quite cool because for scaling laws we can cover a huge range in time, from 10⁻⁵ up to 10.
"

# ╔═╡ 9e714346-0e8d-4de5-85b4-4ede3e275834
"""
	tau_vr(;r₀=171, μ=1/6, γ=γ₀)

Computation of the viscous relaxation time.

#### Mathematics

`` \\tau_{\\text{vr}} = \\frac{r_0\\mu}{\\gamma} ``

where ``r_0`` is the base radius of the droplet, ``\\mu`` is the viscosity and ``\\gamma`` is the surface tension.
"""
function tau_vr(;r₀=171, μ=1/6, γ=γ₀)
	return r₀ * μ / γ
end

# ╔═╡ a2389b72-c460-4026-8676-b3b0fb4acdc2
"""
	tau_ic(;ρ=1, r₀=171, γ=γ₀)

Computation of the inertio-capillary time.

#### Mathematics

`` \\tau_{\\text{ic}} = \\sqrt{\\frac{\\rho r_0^3}{\\gamma}} ``

where ``r_0`` is the base radius of the droplet, ``\\rho`` is the liquids density and ``\\gamma`` is the surface tension.
"""
function tau_ic(;ρ=1, r₀=171, γ=γ₀)
	return sqrt(ρ*r₀^3/γ)
end

# ╔═╡ 2ef048ee-fc8a-4898-a82c-4c33bc1c52b3
md"Loading the data from disk" 

# ╔═╡ ac583c93-ca39-420c-95dc-8db69c55a790
begin
	sys2 = Swalbe.SysConst_1D(L=1024, param=Swalbe.Taumucs(n=9, m=3, Tmax=5000000, δ=5.0))
	df2 = DataFrame()
	# Time
	time = 100:100:size(data)[2]*100
	# Loop over γ
	for i in 1:4
		# Check if there is already a file created
		sim_name = "gamma_$(gamma_labels[i])_tmax_$(sys2.param.Tmax).jld2"
		save_file = string(data_path, sim_name)
		# If so just read it from disc
		if isfile(save_file)
			tmp = load(save_file) |> DataFrame
			tmp.nablaG .= gamma_labels[i]
			df2 = vcat(df2, tmp)
		# If not, compute the evolution of the droplet coalescence
		else
	 		println("Can not find simulation results for run $(save_file)\nCheck the data folder.")
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

However, if we are already taking a look at the data why not ask it more.
We could for example ask if the minimum of the neck is time independent,

$\forall t : h_0(t, x) = h_0(t_0, x_0) ,$

so it stays at the same position throughout the simulation.

Another interesting question to ask is about the symmetry of the neck.
If we were to put a symmetry line along the center of the domain, will the two sides of the neck being mirror symmetric, this is however work in progress.
From the pictures above we already know that this seems to be only the case for the constant surface tension field $\gamma(x) = \gamma_0$.

The is collected into a dataframe which makes it easy to store it and plot specific parts.
"

# ╔═╡ b8eadb42-643c-41e3-ae94-7fe0cefb3d6b
"""
	bridge_height(df; time=100:100:5000000, L=L)

Measurement of bridge height, neck position and skewness.
"""
function bridge_height(df; time=100:100:5000000, L=L, r0=171)
	df_ = DataFrame()
	# Parameter
	center = L÷2
	# Controll values
	time_list = []
	grad_list = []
	# Computed data
	height_list = []
	skew_list = []
	pos_min_list = []
	# Loop through data
	for i in values(gamma_labels)
		tmp = @subset(df, :nablaG .== i) 
		for t in time
			neck_region = tmp[!, Symbol("h_$(t)")][center-r0:center+r0]
			hmin = argmin(neck_region)
			push!(time_list, t)
			push!(grad_list, i)
			push!(height_list, minimum(neck_region))
			push!(pos_min_list, argmin(neck_region))
			push!(skew_list, maximum(reverse(tmp[!, Symbol("h_$(t)")][hmin-100:hmin]) - tmp[!, Symbol("h_$(t)")][hmin+1:hmin+101]))
		end
	end
	df_[!, "nablaG"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_height"] = height_list
	df_[!, "neck_min"] = pos_min_list
	df_[!, "skewness"] = skew_list

	return df_
end

# ╔═╡ 20d4019f-7f79-4f25-bbdd-e439f936fa62
"""
	simple_curvature(x; appr = zeros(length(x)))

Straight forward computation of 

`` \\kappa = \\frac{\\frac{d^2 y}{dx^2}}{[1 + (\\frac{dy}{dx})^2]^{3/2}}``

and returns an array of similar size to the input with curvature values.
"""
function simple_curvature(x; appr = zeros(length(x)))
	p = circshift(x, 1)
    m = circshift(x,-1)
	
    appr .= p .- 2 .* x .+ m
    appr ./= (1 .+ (0.5 .* (p .- m)).^2).^(3/2)

	return appr
end

# ╔═╡ 8adf34e9-7b5c-4b75-8753-2726ffdda706
md"The data is stored in a dataframe with column names:

- **nablaG**: The shape of the surface tension gradient
- **time**: Time steps in lattice Boltzmann units
- **bridge_height**: Minimal thickness between the droplets
- **neck_min**: Position of the minimum on the lattice
- **skewness**: First try to measure asymmetry

Filtering on these columns is rather straight forward with `DataFramesMeta`.
"

# ╔═╡ 2cc159d3-1548-4f46-842d-0ebeecee49be
begin
	analysis_1 = "Neck_bridge_skew.csv"
	frame_analysis_1 = string(data_path, analysis_1)
	if isfile(frame_analysis_1)
		df_secsweep = CSV.File(frame_analysis_1) |> DataFrame
	else
		df_secsweep = bridge_height(df2)
		CSV.write(frame_analysis_1, df_secsweep)
	end
end

# ╔═╡ 6dd016bb-9a59-48ed-93ef-6e5c814037e9
md"Theory suggests that the brdige height should grow according to a powerlaw.
It is therefore necessay to have data for different decades.
With a uniform sampling interval however the density of points looks odd.
Very high densities appear in the late time stages, while data seems to sparse in the early time regime.

There are two ideas to overcome this issue.

1. Do not show points at all
2. Use a log sampling of data points

To normalize the axis for the scaling we compute the normalized time scales and the initial bridge height at $t = 0$.
" 

# ╔═╡ 2e6cd7bf-b650-48d1-aa23-78936404dce1
# Center thickness at t=0 divided by r₀
h0t0 = minimum(h[200:600])/171

# ╔═╡ 9c48f835-3595-4d2d-978f-18fb74b6d3fe
# Time divided by the inertio capillary time
tr = @subset(df_secsweep, :nablaG .== "default").time ./ tau_ic()
#tr_vr = @subset(df_secsweep, :nablaG .== "default").time ./ tau_vr()

# ╔═╡ 2b82f7ec-08de-4cd8-b4ee-7fac37fc4876
"""
	log_t(; tr=tr, log_p_interval=[1, 2, 3, 4, 5, 7])

Data point sampling for log log plots.
"""
function log_t(; tr=tr, log_p_interval=[1, 2, 3, 4, 5, 7])
	some_list = []
	t_len = length(tr)
	j = 1
	while j < t_len + 1
		for i in log_p_interval
			if i * j < t_len
       			push!(some_list, i*j)
			end
    	end
       	j *= 10
    end
	return some_list
end

# ╔═╡ 01c10b58-c4de-456b-a69c-f751eaa879ad
some_list = log_t()

# ╔═╡ 7a9cc2e4-5c43-49e5-abd6-5c614e7b3c30
begin
	p_bridge = plot(tr, @subset(df_secsweep, :nablaG .== "default").bridge_height ./ 171,
	ylabel = "h₀/R₀", 
	xlabel = "t/τ",
	label = "default",
	title = "Evolution bridge height",
	l = (3, :auto),
	xaxis = :log, 
	yaxis = :log,
	xticks=([0.001, 0.01, 0.1, 1, 10], 
	        ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]), # Axis labeling
	legendfontsize = 14,		# legend font size
    tickfontsize = 14,			# tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr, @subset(df_secsweep, :nablaG .== "linear").bridge_height ./ 171, l=(3, :auto),label="linear")
	plot!(tr, @subset(df_secsweep, :nablaG .== "step").bridge_height ./ 171, l=(3, :auto),label="step")
	plot!(tr, @subset(df_secsweep, :nablaG .== "tanh").bridge_height ./ 171, l=(3, :auto),label="tanh")
	fit_t = 2e-4:2e-4:100
	exponent = 2/3
	plot!(fit_t, 0.014 .* fit_t.^(exponent) .+ h0t0, l=(3, :auto, :black), label="∝t^(2/3)")
	plot!(xlim=(5e-4, 30), ylim=(1e-3, 0.1))
end

# ╔═╡ c688523c-2db3-4362-b379-20115cbbfde9
md"Although not perfect but the theory with some fitting seems to capture the trend quite well.
For the fitting we have used

$f(t) = a + b\cdot t^{2/3},$

where 

$a = \frac{h_0(t_0)}{r_0},$ 

and

$b = 0.014.$

Here we hide the data points and instead show a rather smooth curve.
Therefore this is inline with **Do not show data points**

However presentation of data is relatively important.
Doing a simple plot is often not enough to please possible referees.
Because of the fact that the data is discrete it helps to show data points.
In double logarithmic scale this can be a hassle.

That is why in the cell below we try to get a representative subset of data points only, in a log like distrubted fashion.
For that reason we create a time array with appropriate distribution to ensure that the point density does not become unbareable.
Meaning we use **log like display of data points** (there is more data than points!)
"

# ╔═╡ 1c780662-5b9a-4f5d-bb6c-3025bde4e484
md"Using this time subset the data can cleanly displayed as"

# ╔═╡ a1df9654-cafa-4739-be95-d36b0cfa62ac
begin
	logset = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,]
	p_b_simple = plot(tr[some_list], @subset(df_secsweep, :nablaG .== "default").bridge_height[some_list] ./ 171,
	ylabel = "h₀/R₀", 
	xlabel = "t/τ",
	label = "γ(x) = const.",
	# title = "Evolution bridge height",
	m = (9, :auto, 0.6),
	st = :scatter, 				# some recipy stuff
	xaxis = :log, 
	yaxis = :log,
	grid = false,
	xticks=([0.001, 0.01, 0.1, 1, 10], 
	        ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]), # Axis labeling
	legendfontsize = 14,		# legend font size
    tickfontsize = 14,			# tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr[some_list], @subset(df_secsweep, :nablaG .== "step").bridge_height[some_list] ./ 171, 
	m = (9, :auto, 0.6),
	st = :scatter,
	label="γ(x) = Θ(x)"
	)
	plot!(fit_t, 0.0165 .* fit_t.^(exponent), l=(3, :black), label="f(x) ∝ t^(2/3)")
	plot!(xlim=(5e-4, 30), ylim=(1e-3, 0.12))
end

# ╔═╡ 72aba74a-6fd5-4bcc-b0b7-e8b0854b52bd
# Enable to save the figure
# savefig(p_b_simple, "..\\..\\figures\\simple_bridge_height.svg")

# ╔═╡ 32d9caab-4c07-4546-aea7-708fd1fbe6b1
md"The lower plot, while with less data is somewhat better to understand.
The legend here tells us that the blue bullets were created using a constant surface tension throughout the whole domain, while orange squares are generated using a Heaviside function 

$\Theta(x) = \begin{cases}
\gamma_0 \quad\qquad\text{for}\quad x < L/2 \\
\gamma_0 - \epsilon \quad~\text{else}
\end{cases}$

and $f(x)$ here is the same as above with $a = 0$.

The next question we asked the data was if the minimum was moving with time or if its time independent.
The plot below shows what the data has to say about that question.
"

# ╔═╡ ce0bf03e-181b-4abe-89b4-af0ae4217346
begin
	p_pos = plot(tr, @subset(df_secsweep, :nablaG .== "default").neck_min .- 172,
	ylabel = "χ", 
	xlabel = "t/τ",
	label = "default",
	title = "Position of the minimum",
	l=(3, :auto),
	# xaxis=:log, 
	# yaxis=:log,
	legendfontsize = 14,  # legend font size
    tickfontsize = 14,	  # tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr, @subset(df_secsweep, :nablaG .== "linear").neck_min .- 172, l=(3, :auto),label="linear")
	plot!(tr, @subset(df_secsweep, :nablaG .== "step").neck_min .- 172, l=(3, :auto),label="step")
	plot!(tr, @subset(df_secsweep, :nablaG .== "tanh").neck_min .- 172, l=(3, :auto),label="tanh")
end

# ╔═╡ 77ecc143-24c0-448e-93f1-b8cec5b60700
md"One of the interesting things that happen during this experiments is the movement of the minimum.
We define the minimum of the bridge as 

$\chi(t) := pos(h_0) - \frac{L}{2},$

where $pos(g)$ returns the $x$-coordinate of $g(t)$ (somewhat similar to `argmin`).
Clearly for the case of constant γ (blue line) no asymmetry arises and the minimum keeps being pinned in the center thus

$\chi(t)^{\gamma_0} = const.$

Somewhat understandable the minimum moves in the case of a linear surface tension gradient (orange dashed).
The linear surface tension function is with a constant negative slope, therefore the minimum is moving into regions of lower surface tension or

$\chi(t)^{\gamma_l} = kt,$

where $k$ is a positive slope and $t$ is the time.

The smoothed step with a tangent hyperbolicus (violet dashed dotted) yields some interesting dynamics, with an unexpected $\chi(t) = 0$ crossing. 
We will get back to this later.

Almost similar to the constant surface tension, the step function induces a reshaping of the double droplet state and quickly reaches an equilibrium at negativ $\chi$,

$\chi(t)^{\gamma_s} = \chi(t)^{\gamma_0} - l_0,$

where $l_0$ is the shift in distance towards the new static equilibrium.
"

# ╔═╡ 0cf67e97-fb9a-43d2-985f-90e3328f680c
begin
	plot(tr, @subset(df_secsweep, :nablaG .== "default").neck_min .- 172,
	ylabel = "χ", 
	xlabel = "t/τ",
	label = "default",
	title = "Position of the minimum",
	l=(3, :auto),
	st = :samplemarkers, 				# some recipy stuff
	step = 2000, 
	marker = (8, :auto, 0.6),
	legendfontsize = 14,  # legend font size
    tickfontsize = 14,	  # tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr, @subset(df_secsweep, :nablaG .== "linear").neck_min .- 172, l=(3, :auto),label="linear", st = :samplemarkers, marker = (8, :auto, 0.6), step = 2000)
	plot!(tr, @subset(df_secsweep, :nablaG .== "step").neck_min .- 172, l=(3, :auto),label="step", st = :samplemarkers, marker = (8, :auto, 0.6), step = 2000)
	plot!(tr, @subset(df_secsweep, :nablaG .== "tanh").neck_min .- 172, l=(3, :auto),label="tanh", st = :samplemarkers, marker = (8, :auto, 0.6), step = 2000)
end

# ╔═╡ e89f39ed-473c-40aa-ac2e-40844bd80127
md"Pretty much the same with markers for every nth data point.
The discret nature of the plot arises from the fact that we do not interpolate the minimum.
Thus it can only be above a well defined lattice point.

The last point on the first analysis agenda was the skewness.
Does the surface tension gradient induce asymmetry around the center?
The clear answer to this question is, **yes**.
Sadly I have to come up with a better way to measure this.
However below is the first try.
Trying to understand maybe not helpful at all!"

# ╔═╡ 67b0a5a4-f5f3-4778-bc79-49a8e5af9b64
begin
	p_xi = ξ₀ = @subset(df_secsweep, :nablaG .== "default").skewness
	plot(tr, @subset(df_secsweep, :nablaG .== "default").skewness ./ ξ₀,
	ylabel = "ξ/ξ₀", 
	xlabel = "t/τ",
	label = "default",
	title = "Asymmetry around the minimum",
	l=(3, :auto),
	xaxis=:log, 
	# yaxis=:log,
	legendfontsize = 14,  # legend font size
    tickfontsize = 14,	  # tick font and size
    guidefontsize = 15,
	legend=:bottomleft)
	plot!(tr, @subset(df_secsweep, :nablaG .== "linear").skewness ./ ξ₀, l=(3, :auto),label="linear")
	plot!(tr, @subset(df_secsweep, :nablaG .== "step").skewness ./ ξ₀, l=(3, :auto),label="step")
	plot!(tr, @subset(df_secsweep, :nablaG .== "tanh").skewness ./ ξ₀, l=(3, :auto),label="tanh")
	xlims!(1e-2, 50)
end

# ╔═╡ 851750ae-46ec-4853-8c48-06608b2614c5
md"### Results - Playing with the smoothing width

We know from Stefan Krapitschkas experiments that droplets can be kept from coalescing.
In their experiments they use two liquids with different surface tensions, which roughly corresponds to the step surface tension in our simulation.

What happens when the step is relaxed?
Above we saw, the droplets coalesce again.
In fact the only case where the droplets were not coalescing was the step function.
However we like to study the influence of the gradient a little more.

In the following we will only considere a surface tension gradient of typ,

$γ^{smooth}(x) = γ₀[s(x) + (1-s(x))(1-ϵ)],$

with $s(x)$ being

$s(x; a, b) = \Bigg|1 - \left[\frac{1}{2} + \frac{1}{2}\tanh\left(\frac{x - a}{b}\right)\right]\Bigg|,$

and vary the smoothing width $b$.
In the experiments above $b$ was rather large.
We therefore change in the following the value of $b$ because we know in the limit of vanishing width it should more or less be similar to the step case.
" 

# ╔═╡ 5114bb75-637a-4e95-ba3f-e075f181f189
begin
	needed = false
	sys3 = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(n=9, m=3, Tmax=10000000, δ=5.0))
	tanh_label_dict = Dict(1 => "sl_1div80", 2 => "sl_1div90", 3 => "sl_1div100")
	tanh_value_dict = Dict(1 => L÷80, 2 => L÷90, 3 => L÷100)
	df3 = DataFrame()
	for i in 1:3
		ggs = zeros(L)
		gamma_curves_tanh!(ggs, sl=tanh_value_dict[i])
		# Check if there is already a file created
		sim_name_tanh = "gamma_tanh_width_$(tanh_label_dict[i])_tmax_$(sys3.param.Tmax).jld2"
		save_file = string(data_path, sim_name_tanh)
		# If so just read it from disc
		if isfile(save_file)
			tmp = load(save_file) |> DataFrame
			# Add a column for the smoothing width
			tmp.sw .= tanh_value_dict[i]
			df3 = vcat(df3, tmp)
		# If not, compute the evolution of the droplet coalescence
		else
	 		println("Can not find simulation results for run $(save_file)\nCheck the data folder.")
		end
		# Print that you are done
		println("Done with iteration $(tanh_label_dict[i])")
	end
	if needed
		println("Data stored in `df3`.")
	else
		df3 = DataFrame()
		println("For data change `needed` to true at the top of the cell.")
	end
end

# ╔═╡ 10481933-4b7c-4bad-8285-b16bf9c526c3
md"Below you can find data from a run up to 10⁷ time steps, therefore $\approx 45 \frac{t}{\tau}$ in dimensionless time.
Because I am very imperfect human being, I mislabelled the data.
The temporal resolution is not a 100Δt but 200Δt, see the time corrected argument in `bridge_height2`.

The parameter we are most interested in is the bridge height $h_0(t)$.
Let us try to extract the data similar to above using the function `bridge_height2` on `df3`."

# ╔═╡ c6708c4e-1613-44ca-938e-695ae95f7643
md"The data we have collected in `df_tanh` can be used to plot the evolution of the bridge height.
Similar to the plots above, we use a subset of data to have understandable log log plot.
"

# ╔═╡ ba37695d-7f37-4bfe-96b4-118d6253391e
md"Good thing we did not call it a day after the first dataset with the smoothed step surface tension.
In the plot above we have a lot of data.
However, there are two extrema:

1. Constant surface tension: blue bullets
2. Heaviside function: orange squares

The rest of the symbols represent differing smoothing widths.
Interestingly for these rather small smoothing widths the bridge initially starts to grow, but then inverts the trend and seperate to a stable two droplet state.
The meaning of a rather small smoothing width is losely defined as

$b \approx \frac{\max[h(t=0)]}{3}.$

One more observation that is probably not so clear from this data is that saddle point of the curves.
It turns out that there is a correlation between the smoothing width and the saddle point.

Further data is needed to make a more detailed observation of this behaviour.
Luckily there is more data which in the next step will be loaded into a dataframe and analyzed in a similar manner.
"

# ╔═╡ 90d215db-fab4-4b87-9088-6f55e1a8b25d
begin
	sys_long = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(n=9, m=3, Tmax=50000000, δ=5.0, tdump=1000))
	tanh_l_dict = Dict(1 => "sw_20", 2 => "sw_30", 3 => "sw_40", 4 => "sw_50")
	tanh_v_dict = Dict(1 => 20, 2 => 30, 3 => 40, 4 => 50)
	df_l = DataFrame()
	for i in 1:length(tanh_l_dict)  
    	# Check if there is already a file created
    	sim_name_tanh = "gamma_tanh_width_$(tanh_l_dict[i])_tmax_$(sys_long.param.Tmax).jld2"
    	save_file = string(data_path, sim_name_tanh)
		# If so just read it from disc
		if isfile(save_file)
			tmp = load(save_file) |> DataFrame
			# Add a column for the smoothing width
			tmp.sw .= tanh_v_dict[i]
			df_l = vcat(df_l, tmp)
		# If not, compute the evolution of the droplet coalescence
		else
	 		println("Can not find simulation results for run $(save_file)\nCheck the data folder.")
		end
		# Print that you are done
		println("Done with iteration $(tanh_l_dict[i])")
	end
end

# ╔═╡ 56ca66b1-4f03-4132-9693-16a9cb173a24
"""
	bridge_height2(df; time=100:100:5000000, L=L, r0=171)

Measurement of bridge height, neck position and skewness.
Time corrected
"""
function bridge_height2(df; time=100:100:5000000, L=L, r0=171, label_dict=tanh_v_dict)
	df_ = DataFrame()
	# Parameter
	center = L÷2
	# Controll values
	time_list = []
	grad_list = []
	# Computed data
	height_list = []
	skew_list = []
	pos_min_list = []
	# Loop through data
	for i in values(label_dict)
		tmp = @subset(df, :sw .== i) 
		for t in time
			neck_region = tmp[!, Symbol("h_$(t)")][center-r0:center+r0]
			hmin = argmin(neck_region)
			push!(time_list, t * 2)
			push!(grad_list, i)
			push!(height_list, minimum(neck_region))
			push!(pos_min_list, argmin(neck_region))
			push!(skew_list, maximum(reverse(tmp[!, Symbol("h_$(t)")][hmin-100:hmin]) - tmp[!, Symbol("h_$(t)")][hmin+1:hmin+101]))
		end
	end
	df_[!, "width"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_height"] = height_list
	df_[!, "neck_min"] = pos_min_list
	df_[!, "skewness"] = skew_list

	return df_
end

# ╔═╡ bba72789-ce4a-46df-8cc8-3e74476e6dbb
begin
	analysis_tanh = "Neck_bridge_skew_tanh.csv"
	frame_analysis_tanh = string(data_path, analysis_tanh)
	if isfile(frame_analysis_tanh)
		df_tanh = CSV.File(frame_analysis_tanh) |> DataFrame
	else
		df_tanh = bridge_height2(df3, label_dict=tanh_value_dict)
		CSV.write(frame_analysis_tanh, df_tanh)
	end
end

# ╔═╡ f49ce948-6a6f-4418-b8ba-1af6b7c5c695
begin
	tr2 = @subset(df_tanh, :width .== 10).time ./ tau_ic() 
	p_bridge_tanh = plot(tr[some_list], @subset(df_secsweep, :nablaG .== "default").bridge_height[some_list] ./ 171,
	ylabel = "h₀/R₀", 
	xlabel = "t/τ",
	label = "default",
	title = "Evolution bridge height",
	m = (9, :auto, 0.6),
	st = :scatter, 
	xaxis = :log, 
	yaxis = :log,
	xticks=([0.001, 0.01, 0.1, 1, 10], 
	        ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]), # Axis labeling
	legendfontsize = 14,		# legend font size
    tickfontsize = 14,			# tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr[some_list], @subset(df_secsweep, :nablaG .== "step").bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) = Θ(x)")
	plot!(tr2[some_list], @subset(df_tanh, :width .== 10).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=10)")
	plot!(tr2[some_list], @subset(df_tanh, :width .== 11).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=11)")
	plot!(tr2[some_list], @subset(df_tanh, :width .== 12).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=12)")
	# Some fit
	plot!(fit_t, 0.019 .* fit_t.^(exponent), l=(3, :black), label="f(x) ∝ t^(2/3)")
	plot!(xlim=(5e-4, 30), ylim=(1e-3, 0.1))
end

# ╔═╡ 8df94264-e377-4ac3-8866-963c1f97bd35
md"Now that the data is in fetched, lets analyze it and free up the memory."

# ╔═╡ 4e9366a1-c415-4094-b840-0b329f73396c
"""
	bridge_height3(df; time=200:200:10000000, L=L, r0=171)

Measurement of bridge height, neck position and skewness.
"""
function bridge_height3(df; time=200:200:10000000, L=L, r0=171, label_dict=tanh_v_dict)
	df_ = DataFrame()
	# Parameter
	center = L÷2
	# Controll values
	time_list = []
	grad_list = []
	# Computed data
	height_list = []
	skew_list = []
	pos_min_list = []
	# Loop through data
	for i in values(label_dict)
		tmp = @subset(df, :sw .== i) 
		for t in time
			neck_region = tmp[!, Symbol("h_$(t)")][center-r0:center+r0]
			hmin = argmin(neck_region)
			push!(time_list, t)
			push!(grad_list, i)
			push!(height_list, minimum(neck_region))
			push!(pos_min_list, argmin(neck_region))
			# Seems to be something funny happening
			if hmin < 101
				hmin = center
				push!(skew_list, maximum(reverse(tmp[!, Symbol("h_$(t)")][hmin-100:hmin]) - tmp[!, Symbol("h_$(t)")][hmin+1:hmin+101]))
			else
				push!(skew_list, maximum(reverse(tmp[!, Symbol("h_$(t)")][hmin-100:hmin]) - tmp[!, Symbol("h_$(t)")][hmin+1:hmin+101]))
			end
		end
	end
	df_[!, "width"] = grad_list
	df_[!, "time"] = time_list
	df_[!, "bridge_height"] = height_list
	df_[!, "neck_min"] = pos_min_list
	df_[!, "skewness"] = skew_list

	return df_
end

# ╔═╡ 3f3ed980-bd22-455d-9f17-52f589055a82
begin
	ana_tanh = "Neck_bridge_skew_tanh_t5e7.csv"
	frame_tanh_l = string(data_path, ana_tanh)
	if isfile(frame_tanh_l)
		df_tanh_l = CSV.File(frame_tanh_l) |> DataFrame
	else
		df_tanh_l = bridge_height3(df_l, time=1000:1000:50000000, label_dict=tanh_v_dict)
		CSV.write(frame_tanh_l, df_tanh_l)
	end
end

# ╔═╡ 4f2e684d-0420-4edc-9a9f-c7140f05d208
begin
	tr3 = @subset(df_tanh_l, :width .== 20).time ./ tau_ic() 
	p_bridge_tanh_swl = plot(tr[some_list], @subset(df_secsweep, :nablaG .== "default").bridge_height[some_list] ./ 171,
	ylabel = "h₀/R₀", 
	xlabel = "t/τ",
	label = "γ₀",
	title = "Evolution bridge height",
	m = (9, :auto, 0.6),
	st = :scatter, 
	xaxis = :log, 
	yaxis = :log,
	xticks=([0.001, 0.01, 0.1, 1, 10], 
	        ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹"]), # Axis labeling
	legendfontsize = 14,		# legend font size
    tickfontsize = 14,			# tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr[some_list], @subset(df_secsweep, :nablaG .== "step").bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) = Θ(x)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 20).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=20)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 30).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=30)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 40).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=40)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 50).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) ∝ s(x; b=50)")
	# Some fit
	plot!(fit_t, 0.019 .* fit_t.^(exponent), l=(3, :black), label="")
	plot!(xlim=(5e-4, 200), ylim=(1e-3, 0.1))
end

# ╔═╡ a1dfe340-2fe4-498c-828f-3230c4a91a92
md"One more data set to go"

# ╔═╡ 462e0d4c-16f6-4fef-95e2-0cc21fcd5833
begin
	sys_sw = Swalbe.SysConst_1D(L=L, param=Swalbe.Taumucs(n=9, m=3, Tmax=10000000, δ=5.0, tdump=200))
	tanh_l2_dict = Dict(1 => "sw_1", 2 => "sw_2", 3 => "sw_5", 4 => "sw_10")
	tanh_v2_dict = Dict(1 => 1, 2 => 2, 3 => 5, 4 => 10)
	df_sw = DataFrame()
	for i in 1:length(tanh_l2_dict)  
    	# Check if there is already a file created
    	sim_name_tanh = "gamma_tanh_width_$(tanh_l2_dict[i])_tmax_$(sys_sw.param.Tmax).jld2"
    	save_file = string(data_path, sim_name_tanh)
		# If so just read it from disc
		if isfile(save_file)
			tmp = load(save_file) |> DataFrame
			# Add a column for the smoothing width
			tmp.sw .= tanh_v2_dict[i]
			df_sw = vcat(df_sw, tmp)
		# If not, compute the evolution of the droplet coalescence
		else
	 		println("Can not find simulation results for run $(save_file)\nCheck the data folder.")
		end
		# Print that you are done
		println("Done with iteration $(tanh_l2_dict[i])")
	end
end

# ╔═╡ d3bb48e3-550a-4378-a22b-359abb3c494e
begin
	an_tanh = "Neck_bridge_skew_tanh_t1e7.csv"
	frame_tanh_sw = string(data_path, an_tanh)
	if isfile(frame_tanh_sw)
		df_tanh_sw = CSV.File(frame_tanh_sw) |> DataFrame
	else
		df_tanh_sw = bridge_height3(df_sw, time=200:200:10000000, label_dict=tanh_v2_dict)
		CSV.write(frame_tanh_sw, df_tanh_sw)
	end
end

# ╔═╡ 0d58d762-e39b-45a6-907b-105a1d530f76
begin
	tr4 = @subset(df_tanh_sw, :width .== 10).time ./ tau_ic() 
	p_bridge_tanh_sw = plot(tr[some_list], @subset(df_secsweep, :nablaG .== "default").bridge_height[some_list] ./ 171,
	ylabel = "h₀/R₀", 
	xlabel = "t/τ",
	label = "γ₀",
	title = "Evolution bridge height",
	m = (9, :auto, 0.6),
	st = :scatter, 
	xaxis = :log, 
	yaxis = :log,
	xticks=([0.001, 0.01, 0.1, 1, 10, 100], 
	        ["10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹", "10²"]), # Axis labeling
	legendfontsize = 14,		# legend font size
    tickfontsize = 14,			# tick font and size
    guidefontsize = 15,
	legend=:topleft)
	plot!(tr4[some_list], @subset(df_secsweep, :nablaG .== "step").bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="γ(x) = Θ(x)")
	plot!(tr4[some_list], @subset(df_tanh_sw, :width .== 1).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=1)")
	plot!(tr4[some_list], @subset(df_tanh_sw, :width .== 2).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=2)")
	plot!(tr4[some_list], @subset(df_tanh_sw, :width .== 5).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=5)")
	plot!(tr4[some_list], @subset(df_tanh_sw, :width .== 10).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=10)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 20).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=20)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 30).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=30)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 40).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=40)")
	plot!(tr3[some_list], @subset(df_tanh_l, :width .== 50).bridge_height[some_list] ./ 171, m = (9, :auto, 0.6), st = :scatter, label="s(x; b=50)")
	# Some fit
	plot!(fit_t, 0.019 .* fit_t.^(exponent), l=(3, :black), label="")
	plot!(xlim=(5e-4, 200), ylim=(1e-3, 0.1))
end

# ╔═╡ 6ebade34-7de7-4463-a267-5016eb91bedd
md"### Errors

While it may be true that the step at the end of the domain induce some further dynamics, I think it should not appear in this experiment.

See the plot below, having a rather stable double droplet system with one satelite to the right."

# ╔═╡ 47dbb4d1-fada-4c46-acdf-458bee78d9db
begin
	plot(@subset(df_l, :sw .== 20).h_48000000)
	plot!(@subset(df_l, :sw .== 50).h_48000000)
end

# ╔═╡ 29bccf1d-55f8-426a-aab0-59a11eb86a64
"""
	function gamma_curves_tanh_p!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
"""
function gamma_curves_tanh_p!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end

	function smooth_p(l, sl)
		return  (0.5 .+ 0.5 .* tanh.((l .- (L - (2sl+20))) ./ (sl)))
	end
	
	out[:] .= x0 .* smooth(x, L, sl)  .+ (1 .- smooth(x, L, sl)) .* x0 .*(1 - ϵ)
	subspace = L-200:1024
	out[subspace] = x0 .*(1 - ϵ) .* (1 .- smooth_p(subspace, sl)) .+ smooth_p(subspace, sl) .* x0
	return nothing
end

# ╔═╡ fc7b3054-1079-4f37-8526-aee85699a0cc
begin
	ghm = zeros(L) 
	gamma_curves_tanh_p!(ghm, sl=20)
	plot(ghm, legend=:bottomleft)
end

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
# ╠═691876ad-2ee6-4b87-974f-66a3650c4b2f
# ╟─54427765-643f-44fe-84e1-c7c67b2cfe0d
# ╠═eda3dc93-7626-42eb-82a6-b8615bd0f477
# ╠═677ff3cc-4037-4b19-a521-dbca74a635a7
# ╟─c4236e13-fca3-4350-adf7-b98c7bde8a0a
# ╟─3594c0d9-0010-4086-9e7e-163bdf1b0195
# ╟─708f54fc-0bd4-4577-85ec-4faf38029c2f
# ╟─8010c641-a385-4f3f-a88d-817332e45091
# ╟─28025793-b001-4597-aa1c-f2dd06c8a34e
# ╟─76aff625-dcfb-4875-8366-d4c6aac51e54
# ╟─77566ed2-14b6-48c3-831f-80efce8c2e3e
# ╟─ac41b37e-841f-47c6-b5ff-3b10fc2c86ae
# ╟─a0b3c869-3a7e-4b10-a2d8-7021b8c4c54d
# ╟─2f6e1154-eaff-4c20-9fac-290474f45f0b
# ╟─2edc58c6-4ee0-4c5e-8013-311e81820c4c
# ╟─bab2a9ab-1ee1-4146-95e7-03d87a8f9c35
# ╠═35541b43-a834-4f01-ad0d-fa11be9af74b
# ╟─b164e7ec-eaa8-4f2a-b35a-6ac01fd12875
# ╟─9e714346-0e8d-4de5-85b4-4ede3e275834
# ╟─a2389b72-c460-4026-8676-b3b0fb4acdc2
# ╟─2ef048ee-fc8a-4898-a82c-4c33bc1c52b3
# ╟─ac583c93-ca39-420c-95dc-8db69c55a790
# ╟─19772481-02f9-4091-bfc9-e2853e64d4d6
# ╠═0d7c9262-ff25-46bd-9833-bddfd7f958f9
# ╟─50cb438d-501e-411e-844d-e75569ca3c85
# ╟─b8eadb42-643c-41e3-ae94-7fe0cefb3d6b
# ╟─56ca66b1-4f03-4132-9693-16a9cb173a24
# ╟─20d4019f-7f79-4f25-bbdd-e439f936fa62
# ╟─8adf34e9-7b5c-4b75-8753-2726ffdda706
# ╟─2cc159d3-1548-4f46-842d-0ebeecee49be
# ╟─6dd016bb-9a59-48ed-93ef-6e5c814037e9
# ╟─2e6cd7bf-b650-48d1-aa23-78936404dce1
# ╟─9c48f835-3595-4d2d-978f-18fb74b6d3fe
# ╟─2b82f7ec-08de-4cd8-b4ee-7fac37fc4876
# ╟─01c10b58-c4de-456b-a69c-f751eaa879ad
# ╟─7a9cc2e4-5c43-49e5-abd6-5c614e7b3c30
# ╟─c688523c-2db3-4362-b379-20115cbbfde9
# ╟─1c780662-5b9a-4f5d-bb6c-3025bde4e484
# ╟─a1df9654-cafa-4739-be95-d36b0cfa62ac
# ╠═72aba74a-6fd5-4bcc-b0b7-e8b0854b52bd
# ╟─32d9caab-4c07-4546-aea7-708fd1fbe6b1
# ╟─ce0bf03e-181b-4abe-89b4-af0ae4217346
# ╟─77ecc143-24c0-448e-93f1-b8cec5b60700
# ╟─0cf67e97-fb9a-43d2-985f-90e3328f680c
# ╟─e89f39ed-473c-40aa-ac2e-40844bd80127
# ╟─67b0a5a4-f5f3-4778-bc79-49a8e5af9b64
# ╟─851750ae-46ec-4853-8c48-06608b2614c5
# ╠═5114bb75-637a-4e95-ba3f-e075f181f189
# ╟─10481933-4b7c-4bad-8285-b16bf9c526c3
# ╠═bba72789-ce4a-46df-8cc8-3e74476e6dbb
# ╟─c6708c4e-1613-44ca-938e-695ae95f7643
# ╟─f49ce948-6a6f-4418-b8ba-1af6b7c5c695
# ╟─ba37695d-7f37-4bfe-96b4-118d6253391e
# ╟─90d215db-fab4-4b87-9088-6f55e1a8b25d
# ╟─8df94264-e377-4ac3-8866-963c1f97bd35
# ╟─4e9366a1-c415-4094-b840-0b329f73396c
# ╟─3f3ed980-bd22-455d-9f17-52f589055a82
# ╟─4f2e684d-0420-4edc-9a9f-c7140f05d208
# ╟─a1dfe340-2fe4-498c-828f-3230c4a91a92
# ╟─462e0d4c-16f6-4fef-95e2-0cc21fcd5833
# ╟─d3bb48e3-550a-4378-a22b-359abb3c494e
# ╟─0d58d762-e39b-45a6-907b-105a1d530f76
# ╟─6ebade34-7de7-4463-a267-5016eb91bedd
# ╠═47dbb4d1-fada-4c46-acdf-458bee78d9db
# ╠═29bccf1d-55f8-426a-aab0-59a11eb86a64
# ╠═fc7b3054-1079-4f37-8526-aee85699a0cc
# ╟─3091f849-7ce2-400d-9938-4b89215a0bdc
