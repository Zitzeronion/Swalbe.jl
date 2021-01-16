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


function snapshot!(fluid, height, t; dumping = 1000)
    if t % dumping == 0
        fluid[t÷dumping, :] .= vec(Array(height))
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