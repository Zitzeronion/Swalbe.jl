using DrWatson
@quickactivate :Swalbe
using Plots, CUDA, DataFrames, JDF

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
        height, area_lv, area_ls, red_e, h_evo, drop_pos = Swalbe.measure_dropletpatterned(sys, device, radius=85, θₛ=θ_in)
        
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