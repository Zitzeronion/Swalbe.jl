data_path = "..\\..\\data\\Drop_coalescence\\"
L = 1024
x = collect(1:L)
γ = zeros(4,L)
ε = 0.2
γ₀ = 0.0001
γ_bar = (γ₀ + (γ₀ - ε))/2
Δγ = ε
sl = L÷10
gamma_labels = Dict(1 => "default", 2 => "linear", 3 => "step", 4 => "tanh")

"""
	gamma_curves!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)

Different surface tension fields, constant, linear, step and tanh.
"""
function gamma_curves!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	# Smoothing
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	# Constant surface tension
	out[1,:] .= x0
	# Linear depency
	out[2,:] .= x0 .* (1 .- ϵ .* l ./ L)
	# Step function
	out[3,1:L÷2] .= x0
	out[3,L÷2+1:L] .= x0 - x0 * ϵ
	# Tangent hyperbolicus smoothing
	out[4,:] .= x0 .* smooth(x, L, sl) .+ (1 .- smooth(x, L, sl)) .* x0 .*(1 - ϵ) 
	return nothing
end

"""
	gamma_curves_tanh!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)

Different surface tension fields using the tangent hyperbolicus and a varying smoothing width.
"""
function gamma_curves_tanh!(out; x0=γ₀, ϵ=ε, L=L, sl=sl, l=collect(1:length(out)))
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	
	out[:] .= x0 .* smooth(l, L, sl) .+ (1 .- smooth(l, L, sl)) .* x0 .*(1 - ϵ) 
	return nothing
end

"""
	function gamma_curves_tanh_p!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
"""
function gamma_curves_tanh_p!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end

	function smooth_p(l, sl)
		return (0.5 .+ 0.5 .* tanh.((l .- (sl÷2-1)) ./ (sl)))
	end
	
	out[:] .= x0 .* smooth(x, L, sl)  .+ (1 .- smooth_p(x, sl)) .* x0 .*(1 - ϵ) .+ x0 .* smooth_p(x, sl)
	out[L-sl:end] = x0 .*(1 - ϵ) .* smooth_p(x, sl) .+ (1 .- smooth_p(x, sl)) .* x0 .*(1 - ϵ)
	return nothing
end

"""
	plot_data(data; k=1, t1=100, t2=500, t3=1000, leg=true, labels=gamma_labels)

Convienence function to plot the data.
"""
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