"""
    ∇²f!(output, f, γ)
"""
function ∇²f!(output, f, γ)
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
    # In the end it is just a weighted sum...
    output .= -γ .* (2/3 .* (hjp .+ hip .+ him .+ hjm) 
                  .+ 1/6 .* (hipjp .+ himjp .+ himjm .+ hipjm) 
                  .- 10/3 .* f)
    return nothing
end

"""
    ∇f!(outputx, outputy, f)
"""
function ∇f!(outputx, outputy, f)
    # Straight elements j+1, i+1, i-1, j-1
    fip = circshift(f, (1,0))
    fjp = circshift(f, (0,1))
    fim = circshift(f, (-1,0))
    fjm = circshift(f, (0,-1))
    # Diagonal elements  
    fipjp = circshift(f, (1,1))
    fimjp = circshift(f, (-1,1))
    fimjm = circshift(f, (-1,-1))
    fipjm = circshift(f, (1,-1))
    # In the end it is just a weighted sum...
    outputx .= 1/3 .* (fip .- fim) .+ 1/12 .* (fipjp .- fimjp .- fimjm .+ fipjm)
    outputy .= 1/3 .* (fjp .- fjm) .+ 1/12 .* (fipjp .+ fimjp .- fimjm .- fipjm)

    return nothing
end

function ∇f!(outputx, outputy, f, a)
    # Straight elements j+1, i+1, i-1, j-1
    fip = circshift(f, (1,0))
    fjp = circshift(f, (0,1))
    fim = circshift(f, (-1,0))
    fjm = circshift(f, (0,-1))
    # Diagonal elements  
    fipjp = circshift(f, (1,1))
    fimjp = circshift(f, (-1,1))
    fimjm = circshift(f, (-1,-1))
    fipjm = circshift(f, (1,-1))
    # In the end it is just a weighted sum...
    outputx .= a .* (1/3 .* (fip .- fim) .+ 1/12 .* (fipjp .- fimjp .- fimjm .+ fipjm))
    outputy .= a .* (1/3 .* (fjp .- fjm) .+ 1/12 .* (fipjp .+ fimjp .- fimjm .- fipjm))

    return nothing
end