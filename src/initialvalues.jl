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
    singledroplet(T, size(θ,1), size(θ,2), radius, θ, center, device=false)

Generates a fluid configuration of a single droplet in the shape of spherical cap with contact angle `θ`, sphere radius `radius` and centered at `center`.

# Arguments

- `height::Array{undef, 2}`: numerical formate, either `Float64` or `Float32`
- `radius::AbstractFloat`: radius of the undersize(θ,2)ing sphere from which the spherical cap is cut off
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

# ---------------------------------------------------------------- #
# ---------------------- Substrate patterns ---------------------- #
# ---------------------------------------------------------------- #
"""
    boxpattern(θ, θ₀; center=(size(θ,1)÷2, size(θ,2)÷2), δₐ=1/36, side=20)

Defines a quadratic box around a `center` with side length `side`.

# Arguments

- `θ::Array{undef,2}`: Contact angle map
- `θ₀`: Default contact angle `θ`
- `center::Tuple{Int, Int}`: center of the box pattern, default value is `(size(θ,1)÷2, size(θ,2)÷2)`
- `δₐ::Float64`: Contact angle contrast with the substrate, default is set to `1/36` ≈ 5 degrees difference
- `side::Int`: Length of the sides, default is set to `20`

# Examples

```jldoctest
julia> using Swalbe

julia> θ₀ = 1/9;

julia> pattern, polygon = Swalbe.boxpattern(ones(100,100), θ₀);

julia> polygon # Some cool thing we use to create the posize(θ,2)gones, a LazySet
LazySets.VPolygon{Float64,Array{Float64,1}}([[40.0, 40.0], [60.0, 40.0], [60.0, 60.0], [40.0, 60.0]])

julia> pattern[50,50] == pattern[1,1] # In the center there is a different contact angle!
false

```

# References

Realsize(θ,2) not much to say here, check out [LazySets.jl](https://github.com/JuliaReach/LazySets.jl).

"""
function boxpattern(θ, θ₀; center=(size(θ,1)÷2, size(θ,2)÷2), δₐ=1/36, side=20)
    # Vertex vectors
    xcoords = [center[1]-side/2 center[1]-side/2 center[1]+side/2 center[1]+side/2]
    ycoords = [center[2]-side/2 center[2]+side/2 center[2]-side/2 center[2]+side/2]
    # Check if the vertices are on the grid
    for i in 1:length(xcoords)
        if xcoords[i] < 1 || xcoords[i] > size(θ,1)
            xcoords[i] = 1 # If not push them to the origin
        elseif ycoords[i] < 1 || ycoords[i] > size(θ,2)
            ycoords[i] = 1
        end
    end
    # Vertex coordinates
    vertex₁ = [xcoords[1], ycoords[1]]
    vertex₂ = [xcoords[2], ycoords[2]]
    vertex₃ = [xcoords[3], ycoords[3]]
    vertex₄ = [xcoords[4], ycoords[4]]
    # Build a posize(θ,2)gon
    P = VPolygon([vertex₁, vertex₂, vertex₃, vertex₄])
    for i in 1:size(θ,1), j in 1:size(θ,2)
        if [Float64(i), Float64(j)] ∈ P
            θ[j,i] = θ₀ + δₐ
        else
            θ[j,i] = θ₀
        end
    end
    return θ, P
end

"""
    ellipsepattern(θ, θ₀; center=(size(θ,1)÷2, size(θ,2)÷2), δₐ=1/36, a=10, b=5)

Creates an ellipse shaped contact angle defect with contact angle mismatch `δₐ`.

# Arguments

- `θ::Array{undef,2}`: Contact angle map
- `θ₀`: Default contact angle `θ`
- `center::Tuple{Int, Int}`: center of the created pattern, default values is `center = (size(θ,1)÷2, size(θ,2)÷2)`
- `δₐ::Float64`: contact angle mismatch between patch and rest of substrate, default is `δₐ = 1/36` or 5 degrees
- `a::Int`:: semimajor half ax of the ellipse, default value is `a=10`
- `b::Int`:: semiminor half ax of the ellipse, default value is `b=5`

# Examples

```jldoctest
julia> using Swalbe, Test

julia> θ₀ = 1/9;

julia> θ, P = Swalbe.ellipsepattern(ones(100,100), θ₀); # per default the center is in the middle!

julia> @test θ[1,1] == θ₀
Test Passed

julia> @test θ[50,50] == θ₀ + 1/36 # The default increment, is about 5 degrees.
Test Passed

```

# References
Nothing interesting here.

"""
function ellipsepattern(θ, θ₀; center=(size(θ,1)÷2, size(θ,2)÷2), δₐ=1/36, a=10, b=5)
    # Don't know why LazySets needs this, can not just parse center as argument of Ellipsoid constructor
    mid = ones(2)
    mid[1] = center[1]
    mid[2] = center[2]
    # The eigenvalues of the shape matrix, Q in the LazySets docs, needs to be spezified not the axis length
    shape₁ = a^(2)
    shape₂ = b^(2)
    P = Ellipsoid(mid, [shape₁ 0.0; 0.0 shape₂])
    for i=1:size(θ,1), j=1:size(θ,2)
        if [Float64(i), Float64(j)] ∈ P
            θ[j,i] = θ₀ + δₐ
        else
            θ[j,i] = θ₀
        end
    end
    return θ, P
end

"""
    trianglepattern(θ, θ₀; center=(size(θ,1)÷2, size(θ,2)÷2), δₐ=1/36, side=60)

Generates an equilateral triangle centered around `center` with contact angle contrast `δₐ` and side length `side`.

# Arguments
- `θ::Array{undef,2}`: Contact angle map
- `θ₀`: Default contact angle `θ`
- `center::Tuple(Int, Int)`: position of the center of the triangle
- `δₐ::Float64`: contact angle contrast with the rest of the substrate
- `side::Int`: length of the sides of the equilateral triangle

# Examples
```jldoctest
julia> using Swalbe

julia> θ, P = Swalbe.trianglepattern(ones(50,50), 1/9, side=20) # Returns a polygon and the contact angle field
([0.1111111111111111 0.1111111111111111 … 0.1111111111111111 0.1111111111111111; 0.1111111111111111 0.1111111111111111 … 0.1111111111111111 0.1111111111111111; … ; 0.1111111111111111 0.1111111111111111 … 0.1111111111111111 0.1111111111111111; 0.1111111111111111 0.1111111111111111 … 0.1111111111111111 0.1111111111111111], LazySets.VPosize(θ,2)gon{Float64,Array{Float64,1}}([[15.0, 19.226497308103742], [35.0, 19.226497308103742], [25.0, 36.547005383792516]]))

julia> P
LazySets.VPolygon{Float64,Array{Float64,1}}([[15.0, 19.226497308103742], [35.0, 19.226497308103742], [25.0, 36.547005383792516]])

julia> θ[25,25]
0.1388888888888889
```

# References

See also:
"""
function trianglepattern(θ, θ₀; center=(size(θ,1)÷2, size(θ,2)÷2), δₐ=1/36, side=60)
    height = sqrt(3)/2*side
    # Vertices vectors
    xcoords = [center[1]-side/2 center[1]+side/2 center[1]]
    ycoords = [center[2]-height/3 center[2]-height/3 center[2]+2*height/3]
    # Check if the vertices are on the grid
    for i in 1:length(xcoords)
        if xcoords[i] < 1 || xcoords[i] > size(θ,1)
            xcoords[i] = 1 # If not push them to the origin
        elseif ycoords[i] < 1 || ycoords[i] > size(θ,2)
            ycoords[i] = 1
        end
    end
    vertex₁ = [xcoords[1], ycoords[1]]
    vertex₂ = [xcoords[2], ycoords[2]]
    vertex₃ = [xcoords[3], ycoords[3]]
    P = VPolygon([vertex₁, vertex₂, vertex₃])
    for i in 1:size(θ,1), j in 1:size(θ,2)
        if [Float64(i), Float64(j)] ∈ P
            θ[j,i] = θ₀ + δₐ
        else
            θ[j,i] = θ₀
        end
    end
    return θ, P
end

