using DrWatson
@quickactivate :Swalbe
using Plots, CUDA, DataFrames, JDF

# Fluid dynamics we need for the experiment
"""
    run_dropletpatterned()

Simulates an droplet on a patterned substrate
"""
function measure_dropletpatterned(
    sys::Swalbe.SysConst, 
    device::String; 
    radius=20, 
    θ₀=1/6, 
    center=(sys.Lx÷2, sys.Ly÷2), 
    θₛ=ones(sys.Lx, sys.Ly), 
    verbos=true, 
    T=Float64
)
    println("Simulating a droplet on a patterned substrate")
    height_dump = 1000
    df_drop = zeros(sys.Tmax÷height_dump, sys.Lx*sys.Ly) 
    area_lv = zeros(sys.Tmax)
    area_ls = zeros(sys.Tmax)
    drop_dummy = falses(sys.Lx, sys.Ly)
    red_e = zeros(sys.Tmax)
    h_evo = zeros(sys.Tmax)
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
        # Forces are the pressure gradient and the slippage due to substrate liquid boundary conditions
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
        Swalbe.wetted!(area_ls, h_evo, height, t)
        # Far from clean but good enough for now, one should use dedicated arrays in principle.
        Swalbe.surfacearea!(area_lv, red_e, height, θₛ, h∇px, h∇py, dgrad, pressure, t)
        Swalbe.snapshot!(df_drop, height, t, dumping=height_dump)
    end
    return height, area_lv, area_ls, red_e, h_evo, df_drop
end


println("Drop relaxation experiments on patterned substrates")
# Different substrate patches
for shape in ["box" "triangle" "ellipse" "circle"]
    # Different initial volumes
    for R in [80 85 90]
        # Different contact angle mismatch
        for contrast in [1/36 1/18 1/12 1/9 -1/18]
            println("Simulating shape $shape R0 $R contrast $(Int(round(rad2deg(contrast*pi), digits=0)))")
            sys = Swalbe.SysConst(Lx=350, Ly=350, γ=0.01, δ=0.5, n=3, m=2, hmin=0.07, Tmax=1000000)
            df_measures = DataFrame()
            df_drop = DataFrame()
            θₚ = ones(sys.Lx,sys.Ly)
            height_dump = 1000
            boxside = 87
            triside = 132
            ella = 73
            ellb = 33
            rcirc = 49
            # Substrate patterning
            if shape == "box"
                Swalbe.boxpattern(θₚ, 1/9, center=(sys.Lx÷2, sys.Ly÷2), δₐ=-contrast, side=boxside)
            elseif shape == "triangle"
                Swalbe.trianglepattern(θₚ, 1/9, center=(sys.Lx÷2, sys.Ly÷2), δₐ=-contrast, side=triside)
            elseif shape == "ellipse"
                Swalbe.ellipsepattern(θₚ, 1/9, center=(sys.Lx÷2, sys.Ly÷2), δₐ=-contrast, a=ella, b=ellb)
            elseif shape == "ellipse"
                Swalbe.ellipsepattern(θₚ, 1/9, center=(sys.Lx÷2, sys.Ly÷2), δₐ=-contrast, a=rcirc, b=rcirc)
            end
            # Make a cuarray with the substrate pattern
            θ_in = CUDA.adapt(CuArray, θₚ)
            # Actual simulation
            height, area_lv, area_ls, red_e, h_evo, drop_position = measure_dropletpatterned(sys, "GPU", radius=R, θₛ=θ_in)
            println("Writing measurements to dataframes")
            # Filling the dataframes
            df_measures[!, Symbol("A_lv")] = area_lv
            df_measures[!, Symbol("A_sl")] = area_ls
            df_measures[!, Symbol("E_red")] = red_e
            df_measures[!, Symbol("H_max")] = h_evo
            for t in 1:sys.Tmax÷height_dump
                df_drop[!, Symbol("time_step_$(t*height_dump)")] = drop_position[t,:]
            end
            println("Saving dataframe shape $shape R0 $R contrast $(round(contrast, digits=3)) to disk")
            file1 = JDF.save("data/Pattern_substrate/measures_shape_$(shape)_R0_$(R)_contrast_$(Int(round(rad2deg(contrast*pi), digits=0))).jdf", df_measures)
            file2 = JDF.save("data/Pattern_substrate/drop_pos_shape_$(shape)_R0_$(R)_contrast_$(Int(round(rad2deg(contrast*pi), digits=0))).jdf", df_drop)
            # Get back the memory of dataframes
            df_measures = DataFrame()
            df_drop = DataFrame()
            CUDA.reclaim()
            gc()
        end
    end
end
println("Script done, let's have a look at the data :)")

