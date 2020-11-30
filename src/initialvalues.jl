"""
    sinewave1D!(height, h₀, n, ϵ)

Creates a sine wave like height field with `n` full waves, average height `h₀` and displacement magnitude `ϵ`.  
"""
function sinewave1D!(height, h₀, n::Int, ϵ, dims::Int)
    # Some helping dummy array
    hdummy = ones(size(height))
    # Used to perform the accumulate operation
    height .= 1.0
    # This we get n nice sine waves
    hdummy .= 2π .* n .* accumulate(+, height, dims=dims) ./ (size(height,dims)-1)
    # Still need to use the correct undulation ϵ and initial height
    height .= h₀ .* (1.0 .+ ϵ .* (map(sin, hdummy) .- 1/size(height,dims)))
    return nothing
end

"""
    randinterface!(height, h₀, ϵ)

Creates a random height field with average height `h₀` and displacement magnitude `ϵ`.  
"""
function randinterface!(height, h₀, ϵ)
    # Some helping dummy array
    hdummy = ones(size(height))
    # Used to perform the accumulate operation
    height .= 1.0
    # This we get n nice sine waves
    rand!(hdummy) 
    # Still need to use the correct undulation ϵ and initial height
    height .= h₀ .* (1.0 .+ ϵ .* hdummy)
    return nothing
end

"""
    singledroplet(T, lx, ly, radius, θ, center, device=false)

Generates a fluid configuration of a single droplet in the shape of spherical cap with contact angle `θ`, sphere radius `radius` and centered at `center`.

# Arguments

- `height::Array{undef, 2}`: numerical formate, either `Float64` or `Float32`
- `radius::AbstractFloat`: radius of the underlying sphere from which the spherical cap is cut off
- `θ::AbstractFloat`: contact angle in multiples of `π`
- `center::Tuple{Int,Int}`: x and y coordinates of the center of the droplet

# Examples

```jldoctest
julia> using Swalbe, Test

julia> rad = 50; θ = 1/3;

julia> height = Swalbe.singledroplet(ones(100,100), rad, θ, (50,50));

julia> @test maximum(height) == rad * (1 - cospi(θ)) # Simple geometry
Test Passed

julia> argmax(height) # Which is constistent with the center!
CartesianIndex(50, 50)

```

# References

See also: 
"""
function singledroplet(height, radius, θ, center)
    lx, ly = size(height)
    area = 2π * radius^2 * (1- cospi(θ))
    drop = 1 + 0im
    @inbounds for i in 1:lx
        for j in 1:ly
            circ = sqrt((i-center[1])^2 + (j-center[2])^2)
            if circ <= radius
                height[i,j] = (cos(asin(circ/radius)) - cospi(θ)) * radius 
            else
                height[i,j] = 0.05
            end
        end
    end
    @inbounds for i in 1:lx
        for j in 1:ly
            if height[i,j] < 0
                height[i,j] = 0.05
            end
        end
    end
    return height
end