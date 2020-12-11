using DrWatson
using Plots
@quickactivate :Swalbe

sys = Swalbe.SysConst(Lx=100, Ly=100, g=-0.001, γ=0.0005)
h = Swalbe.run_rayleightaylor(sys, "CPU"; h₀=1.0, ϵ=0.01, verbos=true)

heatmap(h[1], aspect_ratio=1, c=:viridis)
gui()
readline()