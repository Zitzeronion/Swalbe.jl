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

Computes a characteristic time scale for an spinodally dewetting film.

# Arguments

- `hᵦ :: Float64`: height at which the disjoining pressure vanishes
- `γ :: Float64`: surface tension value
- `μ :: Float64`: kinematic viscosity, same as dynamic for ρ=1
- `θ :: Float64`: highest contact angle given as radiant, e.g. θ=π/9 for 20 degrees
e
# Mathematics



# Examples
```jldoctest
julia> using Swalbe, Statistics, Test

julia> x = ones(50,50); y = ones(50,50); h = ones(50,50);

julia> Swalbe.thermal!(x, y, h, 0.1, 1/6, 1)

julia> @test mean(x) ≈ 0.0 atol=1e-2
Test Passed

julia> @test mean(y) ≈ 0.0 atol=1e-2
Test Passed

julia> @test var(x) ≈ 2*0.1/11 atol=(2*0.1/11)/10 # var = 2kbt*6*μ/slip
Test Passed

```

# References

- [Mecke, Rauscher](https://iopscience.iop.org/article/10.1088/0953-8984/17/45/042/meta)
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
    surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface,)

Measures the surface area of the liquid vapor interface and the reduced surface energy.

# Arguments

- `area_lv :: Vector{Float64}`: array to store the computed liquid vapor area
- `red_energy :: Vector{Float64}`: array the stores the computed reduced surface energy 
- `height :: Matrix{Float64}`: current height configuration
- `θ :: Matrix{Float64}`: contact angle field distribution
- `∇hx :: Matrix{Float64}`: height gradient x-component 
- `∇hy :: Matrix{Float64}`: height gradient y-component
- `dgrad :: Array{Float64,3}`: dummy array to store derivatives 
- `surface :: Matrix{Float64}`: array that computes locally the liquid vapor surface area 
- `t :: Int`: current time step
- `hthresh :: Float64`: height threshold below which the substrate is considered *dry*

"""
function surfacearea!(area_lv, red_energy, height, θ::Float64, ∇hx, ∇hy, dgrad, surface, t; htresh = 0.055)
    ∇f_simple!(∇hx, ∇hy, height, dgrad)
    surf = 0.0
    surface .= sqrt.(∇hx.^2 .+ ∇hy.^2 .+ 1)
    surf = sum(surface[height .> htresh])
    area_lv[t] = surf
    red_energy[t] = surf - length(height[height .> htresh]) * cospi(θ) 

    return nothing
end

function surfacearea!(area_lv, red_energy, height, θ::Matrix, ∇hx, ∇hy, dgrad, surface, t; htresh = 0.055)
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