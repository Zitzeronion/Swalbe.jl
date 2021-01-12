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
    area_lv = zeros(sys.Tmax)
    area_ls = zeros(sys.Tmax)
    drop_dummy = falses(sys.Lx, sys.Ly)
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
        # Forces are the pressure gradient and the slippage due to substrate liquid boundary conditions
        Fx .= h∇px .+ slipx
        Fy .= h∇py .+ slipy
        Swalbe.equilibrium!(feq, height, velx, vely, vsq)
        Swalbe.BGKandStream!(fout, feq, ftemp, -Fx, -Fy)
        Swalbe.moments!(height, velx, vely, fout)
        Swalbe.wetted!(area_ls, h_evo, height, t)
        Swalbe.fluid_dry!(drop_pos, drop_dummy, height, t)
        # Far from clean but good enough for now, one should use dedicated arrays in principle.
        Swalbe.surfacearea!(area_lv, red_e, height, θₛ, h∇px, h∇py, dgrad, pressure, t)
    end
    return height, area_lv, area_ls, red_e, h_evo, drop_pos
end


println("Starting iteration through surface tensions and sliplengths")
device = "GPU"
θₚ = ones(300,300)
df_measures = DataFrame()
df_drop = DataFrame()
for γ1 in [0.01]
    for δ1 in [0.5]
        println("Simulating with γ: $γ1 and δ: $δ1")
        sys = Swalbe.SysConst(Lx=300, Ly=300, γ=γ1, δ=δ1, n=3, m=2, hmin=0.07, Tmax=10000, tdump=1000)
        Swalbe.boxpattern(θₚ, 1/6, center=(sys.Lx÷2, sys.Ly÷2), δₐ=-1/12, side=40)
        θ_in = CUDA.adapt(CuArray, θₚ)
        height, area_lv, area_ls, red_e, h_evo, drop_pos = measure_dropletpatterned(sys, device, radius=85, θₛ=θ_in)
        
        df_measures[!, Symbol("A_lv")] = area_lv
        df_measures[!, Symbol("A_sl")] = area_ls
        df_measures[!, Symbol("E_red")] = red_e
        df_measures[!, Symbol("H_max")] = h_evo
        df_drop[!, Symbol("pos_drop")] = drop_pos

        CUDA.reclaim()
    end
end

file1 = JDF.save("data/Pattern_substrate/measures.jdf", df_measures)
file2 = JDF.save("data/Pattern_substrate/drop_pos.jdf", df_drop)


