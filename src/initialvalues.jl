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