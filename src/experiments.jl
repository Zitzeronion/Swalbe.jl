"""
    run_dropletpatterned()

Simulates an droplet on a patterned substrate
"""
function measure_dropletpatterned(
    sys::SysConst, 
    device::String; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.Lx÷2, sys.Ly÷2), 
    θₛ=ones(sys.Lx, sys.Ly), 
    verbos=true, 
    T=Float64
)
    println("Simulating a droplet on a patterned substrate")
    area_lv = zeros(sys.Tmax)
    area_ls = zeros(sys.Tmax)
    red_e = zeros(sys.Tmax)
    h_evo = zeros(sys.Tmax)
    drop_pos = falses(sys.Tmax, sys.Lx*sys.Ly)
    fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py = Swalbe.Sys(sys, device, false, T)
    if device == "CPU"
        Swalbe.singledroplet(height, radius, θ₀, center)
    elseif device == "GPU"
        h = zeros(size(height))
        Swalbe.singledroplet(h, radius, θ₀, center)
        height = CUDA.adapt(CuArray, h)
    end
    Swalbe.equilibrium!(fout, height, velx, vely, vsq)
    ftemp .= fout
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        Swalbe.filmpressure!(pressure, height, dgrad, sys.γ, θₛ, sys.n, sys.m, sys.hmin, sys.hcrit)
        Swalbe.∇f!(h∇px, h∇py, pressure, dgrad, height)
        Swalbe.slippage!(slipx, slipy, height, velx, vely, sys.δ, sys.μ)
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
        wetted!(area_ls, drop_pos, h_evo, height, t)
        # surfacearea!(area_lv, red_e, height, θₛ, h∇px, h∇py, dgrad, pressure)
    end
    return height, area_lv, area_ls, red_e, h_evo, drop_pos
end

function wetted!(area_size, drop_pos, maxheight, height, t; hpfilm = 0.055)
    pos = falses(size(height))
    pos[height .> hpfilm] .= true
    area_size[t] = length(findall(height .> hpfilm))
    drop_pos[t,:] .= vec(pos)
    maxheight[t] = maximum(height)

    return nothing
end

function surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface)
    ∇f_simple!(∇hx, ∇hy, height, dgrad)
    surface .=  sqrt.(∇hx.^2 .+ ∇hy.^2 .+ 1)
    surface[height .< 0.055] .= 0.0
    push!(area_lv, sum(surface))
    push!(red_energy, sum(surface) - sum(θ[height .> 0.055]))

    return nothing
end

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