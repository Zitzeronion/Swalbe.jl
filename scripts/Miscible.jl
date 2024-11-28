using Pkg

Pkg.activate(".")

using Swalbe, DelimitedFiles,  DataFrames, Dates,  LaTeXStrings
using Measures, CSV, CUDA, Shell
using Plots; Plots.gr(); ENV["GKSwstype"] = "100"


function run_Miscible(;
    eps=0.01,
    Tmax=Int(1e7),
    dumps=100,
    tdump=Int(floor(Tmax/dumps)),
    L=2^9,
    delta=1.0,
    hmin=0.1,
    hcrit=0.05,
    n=9, m=3,
    gamma= 0.01 .* ones(3,2),
    data=dataall,
    theta=1/9,
    plot=true,
    run=true,
    tau=1.0,
    D=1e-2,
    h_2=3,
    r_1=50, 
    center_weight=0.5
    )
    #Make folders
    datathis= "$(data)/Tmax=$Tmax,L=$L,delta=$delta,gamma=$(replace(replace("$gamma", ";"=>","), " " => "_")),tau=$tau,D=$D,h_2=$h_2,r_1=$r_1,hmin=$hmin,hcrit=$hcrit"
    try
        mkdir(dataall)
	catch
	end
	try
        mkdir(datathis)
	catch
	end
	try
	    mkdir("$(datathis)/figs")
	catch
	end
    try
	    mkdir("$(datathis)/code")
	catch
	end
    #document code used
    cp("scripts/Miscible.jl","$datathis/code/script.jl", force=true)
    hash=Shell.run("git rev-parse HEAD", capture=true)
    open("$datathis/code/cmd.sh", "w") do f
            write(f, "git chekout $hash\n")
	    write(f, "julia script.jl $(ARGS[1]) $(ARGS[2]) \n")
    end
    #Setup system
    sys = Swalbe.SysConstMiscible_1D(L=L, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=delta, n=n, m=m, tau=tau, D=D, hmin=hmin, hcrit=hcrit)
    state = Swalbe.Sys(sys)
    theta_1=acos((sys.gamma[3,1]-sys.gamma[2,1])/sys.gamma[1,1])
    theta_2=acos((sys.gamma[3,2]-sys.gamma[2,2])/sys.gamma[1,2])
    bar_theta=(theta_1+theta_2)/2
    println("theta_1=$theta_1, theta_2=$theta_2, bar_theta=$bar_theta")
    delta_gamma=abs(sys.gamma[1,1]-sys.gamma[1,2])
    bar_gamma=(sys.gamma[1,1]+sys.gamma[1,2])/2
    println("Delta_gamma=$delta_gamma, bar_gamma=$bar_gamma")
    M=3*delta_gamma/(2*bar_gamma*bar_theta^2)
    println("M=$M")
    # state.height[:,1] .= Swalbe.singledroplet(state.height[:,1], r_1, theta, Int(floor(sys.L/2)))
    # state.height[:,2] .= h_2  #.- state.height[:,1] 
    state.height[:,1] .= Swalbe.droplet_base(state.height[:,1], r_1, theta_1/pi, Int(sys.L/2 - r_1/2), precursor=(sys.hmin-2*sys.hcrit)/2)
    state.height[:,2] .= Swalbe.droplet_base(state.height[:,2], r_1, theta_2/pi, Int(sys.L/2 + r_1/2), precursor=(sys.hmin-2*sys.hcrit)/2)
    bridge_height=[]
    i=0
    println("Starting LBM time loop")
    for t in 0:sys.Tmax
        if t % sys.tdump == 0 
            #Store data
	    if run 
            writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_h_t=$t",state.height)
            writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_slip_t=$t",state.slip)
            writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_gamma_t=$t",state.gamma)
            writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_grad_h_t=$t",state.grad_h)
    		else
		    state.height .= readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_h_t=$t")
	    end
	    a, b, mm= try
		    Swalbe.bridge_height(state,sys,shift=12 , prec_fac=10)
	    catch 
		    NaN, NaN, [NaN,NaN]
	    end
            #Plot data
            if plot
		    Plots.plot(state.height[:,1] .+ state.height[:,2], label=L"h",fillalpha=1, fillrange=0, c= RGB.(state.height[:,2]./(state.height[:,2].+state.height[:,1]),0,state.height[:,1]./(state.height[:,2].+state.height[:,1])), title=L"t/t_{dump}=%$i")
                Plots.plot!(state.height[:,1],label="", c=:white, linewidth=4)
                Plots.plot!(state.height[:,1],label=L"h_1", c=:blue, linewidth=2)
                Plots.plot!(state.height[:,2],label="", c=:white,  linewidth=4)
                p1=Plots.plot!(state.height[:,2],label=L"h_2", c=:red,  linewidth=2)
		Plots.plot!([a], seriestype=:vline, label="")
		Plots.plot!([b], seriestype=:vline, label="")
		# Plots.plot!([ra], seriestype=:vline, label="")
		# Plots.plot!([rb], seriestype=:vline, label="")
		Plots.plot!([mm[2]], seriestype=:vline, label="")
		p1=Plots.plot!([mm[1]], seriestype=:hline, label="")
                Plots.plot(state.gamma[:,1], label=L"\gamma")
                Plots.plot!(state.gamma[:,2], label=L"\gamma_2")
                p2=Plots.plot!(state.gamma[:,3], label=L"\gamma_3")
                Plots.plot(p1,p2, plot_title=L"M=%$(round(M,sigdigits=3)), \bar\theta=%$(round(bar_theta, sigdigits=3)), \Delta\gamma=%$(round(delta_gamma, sigdigits=3)), \bar\gamma=%$(round(bar_gamma, sigdigits=3))")
                Plots.savefig("$(datathis)/figs/$(Swalbe.to_4_digits(i)).png")
            end
            #status report
	    bh=mm[1]
	    append!(bridge_height, bh)
	    println("Time step t=$t, mass 1 m_1=$(round(sum(state.height[:,1]),sigdigits=5)), mass 2 m_2=$(round(sum(state.height[:,2]),sigdigits=5)), diff 1 d_1=$(round(maximum(state.height[:,1])-minimum(state.height[:,1]), sigdigits=7)), diff 2 d_2=$(round(maximum(state.height[:,2])-minimum(state.height[:,2]), sigdigits=7)), bridge_height=$(round(bh,sigdigits=3))")
            i+=1
        end
        #actual simulations
        if run 
            # Swalbe.surface_tension!(state,sys)
            Swalbe.surface_tension_smooth!(state,sys, center_weight=center_weight)
            Swalbe.filmpressure!(state, sys)
	    Swalbe.h∇p!(state,sys)
            Swalbe.slippage!(state,sys)
            Swalbe.∇gamma!(state, sys)
	    state.F .= -state.h∇p  .-  state.slip .+ state.stress
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state, sys)
            Swalbe.moments!(state, sys)
        end
    end 
    writedlm("$(datathis)/00000_bridge_height", bridge_height)
    Plots.plot(bridge_height, xlabel=L"t/t_{dump}", ylabel="bridge height")
    # Plots.plot!(t->t, label=L"t", linestyle=:dashdot, color=:grey)
    Plots.savefig("$(datathis)/figs/00000_bridge_height.png")
end



function solve_for_surface_tension(gamma, bar_theta; step=1e-8, sensitivity=0.01, verbose=false)
if verbose 
println("finding surface tensions")
end
	for gamma_sl in gamma[3,1]- gamma[1,1] :step: gamma[3,1]
		gamma[2,1]=gamma_sl
		gamma[2,2]=gamma_sl
		try
		theta_1=acos((gamma[3,1]-gamma[2,1])/gamma[1,1])
    	theta_2=acos((gamma[3,2]-gamma[2,2])/gamma[1,2])
		diff=(theta_1+theta_2)/2-bar_theta
if verbose
		println("gamma_sl=$gamma_sl,theta_1=$theta_1,theta_2=$theta_2,bar_theta=$((theta_1+theta_2)/2),diff=$(diff)")
end
		if abs(diff)<sensitivity
			println("gamma=$gamma")
			return gamma
		end
		catch 
		end
	end
end


# dataall = "data/$(Dates.today())Miscible_coal_params_00"
dataall = "data/$(Dates.today())Miscible_coal_test_00"

gamma=ones(3,2)
gamma_0=0.002
# liquid vapour
# Varyiing contact angles
# DeltaList=LinRange(0.0,1.0,10)
# fixed contact angles
# DeltaList=LinRange(0.0,0.3,10)
DeltaList=LinRange(0.0,0.3,10)
Delta=DeltaList[parse(Int,ARGS[1])]
gamma[1,1]=gamma_0*(1+Delta/2)
gamma[1,2]=gamma_0*(1-Delta/2)
# solid vapour
gamma[3,1]=2*gamma_0
gamma[3,2]=2*gamma_0
# solid liquid
thetaList=LinRange(1/18,1.5/9,10)
theta=thetaList[parse(Int,ARGS[2])]
# theta=1/9
# fixed contact angle
# gamma[2,1]=gamma[3,1]-gamma[1,1]*cospi(theta)
# gamma[2,2]=gamma[3,2]-gamma[1,2]*cospi(theta)
# varying contact angles
solve_for_surface_tension(gamma, theta*pi)
println("gamma=$gamma")
run_Miscible(Tmax=Int(1e7),gamma=gamma, tau=1.0, D=1.5e-4, h_2=5, r_1=200, delta=1.0, hmin=0.2, L=2^10, run=true)
