"""
    filmpressure!(output, f, γ, θ, n, m, hmin, hcrit)
"""
function filmpressure!(output, f, γ, θ, n, m, hmin, hcrit)
    # Straight elements j+1, i+1, i-1, j-1
    hip = circshift(f, (1,0))
    hjp = circshift(f, (0,1))
    him = circshift(f, (-1,0))
    hjm = circshift(f, (0,-1))
    # Diagonal elements  
    hipjp = circshift(f, (1,1))
    himjp = circshift(f, (-1,1))
    himjm = circshift(f, (-1,-1))
    hipjm = circshift(f, (1,-1))
    # Disjoining pressure part
    κ = γ * (1 - cospi(θ)) * (n-1) * (m-1) / ((n-m)*hmin) 

    output .= -γ .* (2/3 .* (hjp .+ hip .+ him .+ hjm) 
                  .+ 1/6 .* (hipjp .+ himjp .+ himjm .+ hipjm) 
                  .- 10/3 .* f) .- κ .* (power_broad.(hmin./(f .+ hcrit), n) .- 
                                         power_broad.(hmin./(f .+ hcrit), m))
    return nothing
    # This computation is correct, at least mathematically!
end

"""
    power_broad(arg, n)

Computes `arg` to the power `n`.


See also: [`filmpressure`](@ref)
"""
function power_broad(arg::Float64, n::Int)
    temp = 1.0
    for i = 1:n
        temp *= arg
    end
    return temp
end

function power_broad(arg::Float32, n::Int)
    temp = 1.0f0
    for i = 1:n
        temp *= arg
    end
    return temp
end