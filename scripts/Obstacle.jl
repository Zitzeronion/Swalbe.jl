using Swalbe, DelimitedFiles

data = "data/24.03.2022_obstacle"


count = "01"
function no_tours_1D()
    try
    mkdir(data)
    catch
    end
    L=2^9
    obst = zeros(L)
    # obst[1]=1
    # obst[2]=1
    # obst[L]=1
    # obst[L-1]=1
    sys = Swalbe.SysConstWithBound_1D{Float64}(obs=obst, L=L,param= Swalbe.Taumucs( Tmax=1000000, n=3, m=2))
    Swalbe.obslist!(sys)
    state = Swalbe.Sys(sys, kind="gamma_bound")
    state.basestate.height .= Swalbe.two_droplets(sys, r₁=144, r₂=144, θ₁=1/9, θ₂=1/9, center=(50,2^9-53)) 
    writedlm("$data/$(count)_initial.csv", state.basestate.height)
    result = Swalbe.time_loop(sys, state, verbose = true)
    writedlm("$data/$(count)_result.csv", result.basestate.height)
    return state, sys
end



no_tours_1D()