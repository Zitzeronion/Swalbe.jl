using Plots, Revise, DataFrames, FileIO, DataFramesMeta, StatsBase, CSV
include("notebooks\\Droplet_coalescence\\helpers.jl")


sys_tanh = Swalbe.SysConst_1D(L=2048, param=Swalbe.Taumucs(n=9, m=3, Tmax=50000000, δ=5.0, tdump=5000))
tanh_l_dict = Dict(1 => "sw_20", 2 => "sw_30", 3 => "sw_40", 4 => "sw_50")
tanh_v_dict = Dict(1 => 20, 2 => 30, 3 => 40, 4 => 50)
tanh_h = 100
data_tanh = zeros(length(tanh_l_dict), sys_tanh.param.Tmax÷sys_tanh.param.tdump, sys_tanh.L)
rad = 500
sim_time_tanh = 5000:5000:50000000
for i in 1:2 # length(tanh_l_dict)
    ggs = zeros(sys_tanh.L)
    # gamma_curves_tanh!(ggs, sl=tanh_v_dict[i])
    if i == 1
        gamma_curves_tanh!(ggs, sl=tanh_h, L=sys_tanh.L)
        sim_name_tanh = "gamma_tanh_width_sw_100_tmax_$(sys_tanh.param.Tmax).jld2"
        # sim_name_tanh = "gamma_tanh_width_$(tanh_l_dict[i])_tmax_$(sys_tanh.param.Tmax).jld2"
    else
        ggs .= 0.0001
        sim_name_tanh = "gamma_tanh_width_sw_0_tmax_$(sys_tanh.param.Tmax).jld2"
    end
    # Check if there is already a file created
    save_file = string(data_path, sim_name_tanh)
    data_tanh[i, :, :] = Swalbe.run_gamma(sys_tanh, ggs, r₁=rad, r₂=rad, dump=sys_tanh.param.tdump, verbos=false)
    # Think about storing them on the disc
    df_fluid = Dict()
    # Loop through the time once more
    for t in 1:size(data_tanh)[2]
        df_fluid["h_$(sim_time_tanh[t])"] = data_tanh[i,t,:]
    end
    save(save_file, df_fluid)
end