"""
    wetted!(area_size, drop_pos, maxheight, height, t; hthresh = 0.055)

Measures the wetted area and maximal height of the film at time step `t`. 
"""
function wetted!(area_size, maxheight, height, t; hthresh = 0.055)
    area_size[t] = length(findall(height .> hthresh))
    maxheight[t] = maximum(height)

    return nothing
end

"""
    fluid_dry!(fluid, dummy, height, t; hthresh = 0.055)

Tracks the location of the thin film as a boolean field.
"""
function fluid_dry!(fluid, dummy, height, t; hthresh = 0.055)
    dummy .= false
    dummy[height .> hthresh] .= true
    fluid[t, :] .= vec(dummy)
    
    return nothing
end

"""
    t0(;hᵦ=0.07, γ=0.01, μ=1/6, θ=1/6)

Computes a characteristic time scale for an unstable film.
"""
function t0(;hᵦ=0.07, γ=0.01, μ=1/6, θ=1/6)
    qsq = hᵦ * (1 - cospi(θ)) * (2 - 3 * hᵦ) 
    charT = 3 * μ / (γ * qsq^2)

    return charT
end

"""
    snapshot!(snap, field, t; dumping)

Makes a copy of the current state of an input array `in` and saves the vectorized values as a column to `out`.

Function that fills a preallocated array out with a time series of system snap shots of e.g. the height field `h`.

# Arguments

- `snap :: Array{Number,2}`: Array that stores the snap shots as columns
- `field :: Array{Number,2}`: Input argument, e.g. ``h(\\mathbf{x},t)``
- `t :: Int`: The current time step
- `dumping :: Int`: Sampling frequency

# Examples

```jldoctest
julia> using Swalbe, Test

julia> h1 = reshape(collect(1:25),5,5); h2 = reshape(collect(5:5:125),5,5);

julia> snapshot = zeros(2, 25);

julia> Swalbe.snapshot!(snapshot,h1,10,dumping=10)

julia> Swalbe.snapshot!(snapshot,h2,20,dumping=10)

julia> @test all(h1 .== reshape(snapshot[1,:],5,5))
Test Passed

julia> @test all(h2 .== reshape(snapshot[2,:],5,5))
Test Passed

```
# References

See also: The [scripts folder](https://github.com/Zitzeronion/Swalbe.jl/tree/master/scripts) 
"""
function snapshot!(snap, field, t; dumping = 1000)
    if t % dumping == 0
        snap[t÷dumping, :] .= vec(Array(field))
    end
    
    return nothing
end

"""
    surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface)

Measures the surface area of the liquid vapor interface and the reduced surface energy.
"""
function surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface, t; htresh = 0.055)
    ∇f_simple!(∇hx, ∇hy, height, dgrad)
    surf = 0.0
    surface .= sqrt.(∇hx.^2 .+ ∇hy.^2 .+ 1)
    surf = sum(surface[height .> htresh])
    area_lv[t] = surf
    red_energy[t] = surf - sum(cospi.(θ[height .> htresh]))

    return nothing
end

"""
    ∇f_simple!(outputx, outputy, f, dgrad) 

Simple gradient calculation for the differential surface area.
"""
function ∇f_simple!(outputx, outputy, f, dgrad)
    fip, fjp, fim, fjm, fipjp, fimjp, fimjm, fipjm = Swalbe.viewneighbors(dgrad)
    # Straight elements  j+1, i+1, i-1, j-1
    circshift!(fip, f, (1,0))
    circshift!(fjp, f, (0,1))
    circshift!(fim, f, (-1,0))
    circshift!(fjm, f, (0,-1))
    # Diagonal elements  
    circshift!(fipjp, f, (1,1))
    circshift!(fimjp, f, (-1,1))
    circshift!(fimjm, f, (-1,-1))
    circshift!(fipjm, f, (1,-1))
    # In the end it is just a weighted sum...
    outputx .= -1/3 .* (fip .- fim) .- 1/12 .* (fipjp .- fimjp .- fimjm .+ fipjm)
    outputy .= -1/3 .* (fjp .- fjm) .- 1/12 .* (fipjp .+ fimjp .- fimjm .- fipjm)

    return nothing
end