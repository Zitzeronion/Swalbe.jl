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
    surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface)

Measures the surface area of the liquid vapor interface and the reduced surface energy.
"""
function surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface)
    ∇f_simple!(∇hx, ∇hy, height, dgrad)
    surface .=  sqrt.(∇hx.^2 .+ ∇hy.^2 .+ 1)
    surface[height .< 0.055] .= 0.0
    push!(area_lv, sum(surface))
    push!(red_energy, sum(surface) - sum(θ[height .> 0.055]))

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