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
This behavior can be explained following a self-similar ansatz for the thickness of the touching region, the socalled bridge height $h_0$.
In the low contact angle regime one finds the solution $h_0(t) \propto t^{2/3}$.
We gona take this case as starting point and introduce a surface tension gradient to the problem.
This gradient can for example be introduce with the addition of surfactants to only one of the two droplets.

## Numerical setup

Simulations or numerical experiments are performed using the solver called *[Swalbe.jl](https://github.com/Zitzeronion/Swalbe.jl)*.
As of now this simulations are work in progress and surface tension gardients are not yet part of master, therefore we use *[DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl)* to  call the local version of *Swalbe*.
For data analysis purposes two more packages are added.
On the one hand *[Plots.jl](https://github.com/JuliaPlots/Plots.jl)* for basic plotting features and *[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)* for data storage and analysis.
First *DrWatson* is loaded to access a local version of *Swalbe* and subsequently we load *Plots*,*DataFrames* and *FileIO* to write data to the disk."

# ╔═╡ 54427765-643f-44fe-84e1-c7c67b2cfe0d
md"## Surface tension

A lot of research has been dedicated towards the understanding of droplet coalescence.
However as far the authors know the case with a variable surface tension along the coalescence axis has not been studied yet.
This case can arise for example assuming that one of the droplets has a non vanishing surfactant concentration or one uses a well defined heating source to locally heat the substrate (assuming that temperature distribution is fast as compared to the coalescence dynamics).
Both scenarios nevertheless introduce a surface tension gardient.

We are going to study the effect of this gradient using a simple one dimensional model.
The concrete choice for the spatially resolved surface tension is given by either one of three functions,

$γ^{lin.}(x) = γ₀\Bigg(1 - ϵ \frac{x}{L}\Bigg),$

where γ₀ is a characteristic surface tension value, ϵ is a small deviation and $L$ is the size along the system. 
Besides this linear model another surface tension function is given as

$γ^{heavi}(x) = \begin{cases}γ₀\qquad\qquad\text{for x < L/2}\\γ₀(1-ϵ)\quad \text{else}\end{cases},$

with a Heaviside step function at the touching point.
A third option is given using a tagent hyperbolicus to smoothen out the discontinuity of the Heaviside function. 
We use 

$γ^{smooth}(x) = γ₀[s(x) + (1-s(x))(1-ϵ)],$

with 

$s(x) = \Bigg|1 - \left[\frac{1}{2} + \frac{1}{2}\tanh\left(\frac{x - a}{b}\right)\right]\Bigg|$
"

# ╔═╡ eda3dc93-7626-42eb-82a6-b8615bd0f477
begin
	L = 1024
	x = zeros(L)
	γ = zeros(3,L)
	smooth = abs.(1 .- (0.5 .+ 0.5 .* tanh.((collect(1:L) .- L÷2) ./ (L÷10))))
end

# ╔═╡ bcb187be-9a59-46ce-aebe-74e7003077d8
function gamma_curves!(x; x0=0.01, ϵ=0.1)
	x[1,:] .= x0 .* (1 .- ϵ .* collect(1:L) ./ L)
	x[2,1:L÷2] .= x0
	x[2,L÷2+1:L] .= x0 - x0 * ϵ
	# x[3,:] .= x0 .* exp.(-collect(1:L)/L)
	x[3,:] .= x0 .* smooth .+ (1 .- smooth) .* x0 .*(1 - ϵ) 
end

# ╔═╡ 677ff3cc-4037-4b19-a521-dbca74a635a7
gamma_curves!(γ)

# ╔═╡ 708f54fc-0bd4-4577-85ec-4faf38029c2f
begin
	plot(collect(1:L), γ[1,:], w=3, label="lin", xlabel="x/Δx", ylabel="γ")
	plot!(collect(1:L), γ[2,:], w=3, label="step")
	plot!(collect(1:L), γ[3,:], w=3, label="tanh")
end

# ╔═╡ 7b53009d-9000-4885-8fdf-540eb93bc7d9
begin
	daty = Swalbe.SysConst_1D(L=1024)
	daty.θ * 180
end

# ╔═╡ 09a80dac-0cd5-42f3-9676-e412a58f58db
	begin
	# The initial configuration of the numerical experiment
	rad = 500
	h = Swalbe.two_droplets(Swalbe.SysConst_1D(L=1024), r₁=rad, r₂=rad)
	x_ = 1:1024
	p0 = plot(x_ , h, w=3, aspect_ratio=7, label="Initial conf.", xlabel="x [Δx]", ylabel="h(x)")
	ylims!(0,40)
end

# ╔═╡ Cell order:
# ╟─bb534270-0e59-4c41-a825-fd6dc0fb4a7e
# ╠═52725098-857a-4301-b86b-d9cd819de541
# ╠═005b7402-df13-410c-88e1-69fec65b5ff0
# ╠═f75e02be-8f28-11ec-16ae-8bc62cd08643
# ╠═54427765-643f-44fe-84e1-c7c67b2cfe0d
# ╠═eda3dc93-7626-42eb-82a6-b8615bd0f477
# ╠═bcb187be-9a59-46ce-aebe-74e7003077d8
# ╠═677ff3cc-4037-4b19-a521-dbca74a635a7
# ╠═708f54fc-0bd4-4577-85ec-4faf38029c2f
# ╠═7b53009d-9000-4885-8fdf-540eb93bc7d9
# ╠═09a80dac-0cd5-42f3-9676-e412a58f58db
