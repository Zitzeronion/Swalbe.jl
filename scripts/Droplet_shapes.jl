using DrWatson
@quickactivate :Swalbe
using Plots, CUDA, CSV, DataFrames

println("Starting iteration through surface tensions and sliplengths")
df_area = DataFrame()
df_height = DataFrame()
for γ1 in [0.005]
    for δ1 in [1.0 0.5 0.25]
        println("Simulating with γ: $γ1 and δ: $δ1")
        itname = Symbol("ga_$(γ1)_sl_$(δ1)")
        sys = Swalbe.SysConst(Lx=300, Ly=300, γ=γ1, δ=δ1, n=3, m=2, hmin=0.07, Tmax=100000, tdump=1000000)
        h, diff = Swalbe.run_dropletrelax(sys, "GPU", radius=85)
        # rads = sqrt.(diff ./ π)
        # p = plot(rads/rads[1], xaxis=:log, xlabel="time", ylabel="R/R₀")
        # savefig(p, "../images/rads_$(γ1)_$(δ1).pdf")
        # pp = heatmap(Array(h), aspect_ratio=1, c=:viridis)
        # savefig(pp, "../images/shape_$(γ1)_$(δ1).png")
        df_area[!, itname] = diff
        df_height[!, itname] = vec(Array(h))
        # @tagsave(datadir("Relax_droplet", savename(itname, "bson")),
        #          @dict diff vec(h) itname)
        CUDA.reclaim()
    end
end

CSV.write("data/Relax_droplet/Area_evolution_smaller_g.csv", df_area)
CSV.write("data/Relax_droplet/Height_final_smaller_g.csv", df_height)

# ;gzip data/Relax_droplet/Area_evolution.csv
# ;gzip data/Relax_droplet/Height_final.csv
