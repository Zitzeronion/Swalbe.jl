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
function gamma_curves_tanh!(out; x0=γ₀, ϵ=ε, l=x, L=L, sl=sl)
	function smooth(l, L, sl)
		return abs.(1 .- (0.5 .+ 0.5 .* tanh.((l .- L÷2) ./ (sl))))
	end
	
	out[:] .= x0 .* smooth(x, L, sl) .+ (1 .- smooth(x, L, sl)) .* x0 .*(1 - ϵ) 
	return nothing
end