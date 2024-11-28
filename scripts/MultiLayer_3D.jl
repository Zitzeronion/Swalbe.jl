using Pkg

Pkg.activate(".")

using Swalbe, DelimitedFiles,  DataFrames, Dates,  LaTeXStrings
using Measures, CSV, CUDA, Shell, FFTW, LinearAlgebra, Statistics #, Images, ImageFiltering
using Plots; Plots.gr(); ENV["GKSwstype"] = "100"


function gamma_to_string(gamma)
    return "$(round(gamma[1,2],sigdigits=3))_$(round(gamma[1,3],sigdigits=3))_$(round(gamma[1,4],sigdigits=3))_$(round(gamma[2,3],sigdigits=3))_$(round(gamma[2,4],sigdigits=3))_$(round(gamma[3,4],sigdigits=3)),spread_coef=$(round(gamma[1,3]-gamma[1,2]-gamma[2,3],sigdigits=3))_$(round(gamma[2,4]-gamma[2,3]-gamma[3,4],sigdigits=3))_$(round(gamma[2,3]+gamma[1,4]-gamma[2,4]-gamma[1,3], sigdigits=3))"
end

function slip_coef_to_string(slip_coef)
    return "$(round(slip_coef[1,1,], sigdigits=3))_$(round(slip_coef[1,2], sigdigits=3))_$(round(slip_coef[2,1], sigdigits=3))_$(round(slip_coef[2,2],sigdigits=3))"
end

function mu_to_string(mu)
    return "$(round(mu[1],sigdigits=3))_$(round(mu[2], sigdigits=3))"
end






function run_multilayer(;
    h_1=2,
    h_2=1,
    eps=[0.001,0.001],
    Tmax=Int(1e7),
    dumps=100,
    tdump=Int(floor(Tmax/dumps)),
    Lx=2^8,
    Ly=2^8,
    delta=1.0,
    delta_2=1.0,
    n=9, m=3,
    gamma= [0 0.001 0.001 0.001 ; 0 0 0.001 0.001 ; 0 0 0 0.001; 0 0 0 0],
    data=dataall,
    plot=true,
    run=true,
    hmin=0.1,
    mu=[1/6 1/6],
    device="GPU", 
    devicenumber=1,
	uselogo=false
    )
    CUDA.device!(devicenumber)
    #Setup system
    sys = Swalbe.SysConstMultiLayer(Lx=Lx, Ly=Ly, Tmax=Tmax, tdump=tdump, gamma=gamma, delta=[delta,delta_2], n=n, m=m, hmin=hmin, mu=mu)
    state = Swalbe.Sys(sys, device=device)
    #Make folders
    datathis= "$(data)/h_1=$h_1,h_2=$h_2,Tmax=$Tmax,Lx=$Lx,Ly=$Ly,delta=$(delta)_$(delta_2),gamma=$(gamma_to_string(gamma)),n=$n,m=$m,hmin=$hmin,mu=$(mu_to_string(sys.mu)),eps=$(eps[1])_$(eps[2])"
    host_height=zeros(Lx,Ly,2)
    host_p=zeros(Lx,Ly,2)
    host_grad_h_sqd=zeros(Lx,Ly,1)
    host_Fx=zeros(Lx,Ly,2)
    host_Fy=zeros(Lx,Ly,2)
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
    println("$((sys.gamma[1,3]-sys.gamma[1,2]-sys.gamma[2,3]))")
    println("$(sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4])")
    println("$(sys.gamma[2,3]+sys.gamma[1,4]-sys.gamma[2,4]-sys.gamma[1,3])")
    #document code used
    cp("scripts/MultiLayer_3D.jl","$datathis/code/script.jl", force=true)
    hash=Shell.run("git rev-parse HEAD", capture=true)
    open("$datathis/code/cmd.sh", "w") do f
            write(f, "git chekout $hash\n")
            write(f, "julia script.jl  \n")
    end
    #initial conditions
    if Ly>1
    	host_height[:,:,1] .= Swalbe.browniannoise2D!(host_height[:,:,1], h_1, eps[1], Lx/2-Lx/10)
    	host_height[:,:,2] .= Swalbe.browniannoise2D!(host_height[:,:,2], h_2, eps[2], Lx/2-Lx/10)
    else
    	host_height[:,1,1] .= Swalbe.browniannoise!(host_height[:,1,1], h_1, eps[1], Lx/2-Lx/10)
    	host_height[:,1,2] .= Swalbe.browniannoise!(host_height[:,1,2], h_1, eps[1], Lx/2-Lx/10)
    end
    # host_height[:,:,1] .= [h_1 + 0.25*sin(i*pi*2/Lx) for i in 1:Lx, j in 1:Ly]
    # # host_height[:,:,2] .= [h_2 + 0.25*sin(i*pi*2/Lx) for i in 1:Lx, j in 1:Ly]
    # host_height[:,:,2] .= h_2
    CUDA.copyto!(state.height, host_height)
    # state.height[:,1] .= h_1
    i=0
    println("Starting LBM time loop")
    for t in 0:sys.Tmax
        if t % sys.tdump == 0 
            #Store data
            if run 
                CUDA.copyto!(host_height,state.height)
                CUDA.copyto!(host_p,state.pressure)
                CUDA.copyto!(host_Fx,state.Fx)
                CUDA.copyto!(host_Fy,state.Fy)
                CUDA.copyto!(host_grad_h_sqd,state.grad_h_sq)
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_h1_t=$t",host_height[:,:,1])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_h2_t=$t",host_height[:,:,2])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_p1_t=$t",host_p[:,:,1])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_p2_t=$t",host_p[:,:,2])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fx1_t=$t",host_Fx[:,:,1])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fx2_t=$t",host_Fx[:,:,2])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fy1_t=$t",host_Fy[:,:,1])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fy2_t=$t",host_Fy[:,:,2])
                writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_grad_h_sq_t=$t",host_grad_h_sqd[:,:,1])
            end
            # writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_slipx_t=$t",state.slipx)
            # writedlm("$(datathis)/$(Swalbe.to_4_digits(i))_slipy_t=$t",state.slipy)
            #Plot data
            if plot
                host_height[:,:,1]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_h1_t=$t")
                host_height[:,:,2]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_h2_t=$t")
                host_p[:,:,1]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_p1_t=$t")
                host_p[:,:,2]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_p2_t=$t")
                host_Fx[:,:,1]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fx1_t=$t")
                host_Fx[:,:,2]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fx2_t=$t")
                host_Fy[:,:,1]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fy1_t=$t")
                host_Fy[:,:,2]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_Fy2_t=$t")
                host_grad_h_sqd[:,:,1]=readdlm("$(datathis)/$(Swalbe.to_4_digits(i))_grad_h_sq_t=$t")
                clims=(minimum(host_height[:,:,1]), maximum(host_height[:,:,1] .+ host_height[:,:,2]))
                p1=Plots.heatmap(host_height[:,:,1], title=L"z_1", c=:darktest, clims=clims)
                p2=Plots.heatmap(host_height[:,:,1] .+ host_height[:,:,2], title=L"z_2", c=:lighttest, clims=clims)
                p1_1=Plots.surface(host_height[:,:,1], c=:darktest, clims=clims)
                p1_1=Plots.surface!(host_height[:,:,1] .+ host_height[:,:,2], c=:lighttest, camera=(50,40), clims=clims, cbar=false)
                Plots.plot(host_height[:,Int(sys.Ly/2),1] .+ host_height[:,Int(sys.Ly/2),2], label="",fillalpha=1, fillrange=0, c=:cyan)
                p3=Plots.plot!(host_height[:,Int(sys.Ly/2),1], label="",fillalpha=1, fillrange=0, c=:steelblue)
                Plots.plot(host_p[:,1,1], label=L"p_1")
                p4=Plots.plot!(host_p[:,1,2], label=L"p_2")
                p5=Plots.plot(host_grad_h_sqd[:,1,1], label=L"(\nabla h )^2")
                Plots.plot(host_Fx[:,1,1], label="Fx1")
                Plots.plot!(host_Fx[:,1,2], label="Fx2")
                Plots.plot!(host_Fy[:,1,2], label="Fy2")
                p6=Plots.plot!(host_Fy[:,1,1], label="Fy1")
                Plots.plot(p1_1, p3, p1, p2, size=(1000,800), plot_title=L"t/t_{dump}=%$i")
                Plots.savefig("$(datathis)/figs/$(Swalbe.to_4_digits(i)).png")
            end
            #status report
            println("Time step t=$t, mass 1 m_1=$(round(sum(state.height[:,:,1]),sigdigits=5)), mass 2 m_2=$(round(sum(state.height[:,:,2]),sigdigits=5)), diff 1 d_1=$(round(maximum(state.height[:,:,1])-minimum(state.height[:,:,1]), sigdigits=7)), diff 2 d_2=$(round(maximum(state.height[:,:,2])-minimum(state.height[:,:,2]), sigdigits=7))")
            i+=1
        end
        #actual simulations
        if run 
            Swalbe.filmpressure!(state, sys)
	    Swalbe.h∇p!(state,sys)
            Swalbe.slippage!(state, sys)
            # Swalbe.slippage_no_slip!(state, sys)
            # Swalbe.slippage_fake!(state,sys)
	    state.Fx .= -state.h∇px  .+  state.slipx
	    state.Fy .= -state.h∇py  .+  state.slipy
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state)
            Swalbe.moments!(state, sys)
        end
    end 
end


# In order to fulfill the LSA mu cannot be samller then 0.1, the single layer has similar restriction just not that strict, 
# I'd argue that this is due to the leaving of the Stokes regime. Whether the Numerical results are valid outside the Stokes 
# regime is an open question. I'd argue yes but that is hard to check. 

# dataall = "data/$(Dates.today())MultiLayer_sin_wave_test_01"
# dataall = "data/$(Dates.today())MultiLayer_GPU_first_showcase"



# Lower layer doesn't break
dataall = "data/$(Dates.today())_3D_test_00"
gamma_0=0.01
gamma=gamma_0*ones(4,4)
gamma[2,4]=gamma_0*(1+cospi(2/9))
gamma[1,3]=gamma_0*(2.1)
gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]
run_multilayer(Tmax=Int(1e7),gamma=gamma, h_2=1.0, h_1=6.0, delta=1.0, delta_2=1.0, Lx=2^8,Ly=2^8, mu=[0.1 0.1], hmin=0.2, devicenumber=0)


# # Upper layer doesn't break
# gamma_0=0.01
# gamma=gamma_0*ones(4,4)
# gamma[2,4]=gamma_0*(2.1)
# gamma[1,3]=gamma_0*(1+cospi(2/9))
# gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]
# run_multilayer(Tmax=Int(1e6),gamma=gamma, h_2=1, h_1=1, delta=2.0, delta_2=2.0, Lx=2^8, Ly=2^8, mu=[0.1 0.1], hmin=0.2, devicenumber=3)


# # Both don't break
# gamma_0=0.01
# gamma=gamma_0*ones(4,4)
# gamma[2,4]=gamma_0*(2)
# gamma[1,3]=gamma_0*(2)
# gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]
# run_multilayer(Tmax=Int(1e4),gamma=gamma, h_2=1, h_1=1, delta=2.0, delta_2=2.0, Lx=2^8, Ly=2^8, mu=[0.1 0.1], hmin=0.2)

# # Drop covers drop
# gamma_0=0.01
# gamma=gamma_0*ones(4,4)
# gamma[2,4]=gamma_0*(1+cospi(1/18))
# gamma[1,3]=gamma_0*(1+cospi(2/9))
# gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3] -0.05*gamma_0
# run_multilayer(Tmax=Int(1e7),gamma=gamma, h_2=2, h_1=1.0, delta=2.0, delta_2=2.0, L=2^8, mu=[0.1 0.1], hmin=0.2)

# # drop on drop
# gamma_0=0.001
# gamma=gamma_0*ones(4,4)
# gamma[2,3]=gamma_0*(2)
# gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]
# run_multilayer(Tmax=Int(1e6),gamma=gamma, h_2=2, h_1=2, delta=2.0, delta_2=2.0, Lx=2^8, Ly=2^8, mu=[0.1 0.1], hmin=0.2, run=true, plot=true, device="GPU", devicenumber=1)




# # drop besides drop
# # First try
# gamma_0=0.02
# gamma=gamma_0 .* ones(4,4)
# gamma[2,4]=gamma_0*(1+cospi(2/9))
# gamma[1,3]=gamma_0*(1+cospi(2/9))
# gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]+0.1*gamma_0
# # Second try
# # gamma=gamma_0*ones(4,4)
# # gamma[2,3]=gamma_0*(2)
# # gamma[1,3]=gamma_0*(1+cospi(1/9))
# # gamma[1,2]=gamma_0*0.1*cospi(1/9)
# # gamma[2,4]=gamma_0*(1+cospi(1/9))
# # gamma[3,4]=gamma_0*0.1*cospi(1/9)
# # gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]
# run_multilayer(Tmax=Int(1e7),gamma=gamma, h_2=2, h_1=2, delta=2.0, delta_2=2.0, Lx=2^8, Ly=2^8, mu=[0.1 0.1], hmin=0.2, devicenumber=2)


# # Sometthing else
# th_list=[0,1/18,1/9,2/9]
# e_list=[1.0,0.1,0,-0.1,-1.0]
# theta_1=th_list[parse(Int,ARGS[1])]
# theta_2=th_list[parse(Int,ARGS[2])]
# gamma_0=0.01
# extra=e_list[parse(Int,ARGS[3])]
# gamma=gamma_0*ones(4,4)
# gamma[2,4]=gamma_0*(1+cospi(theta_2))
# gamma[1,3]=gamma_0*(1+cospi(theta_1))
# gamma[1,4]=gamma[2,4]+gamma[1,3]-gamma[2,3]+gamma_0*extra
# run_multilayer(Tmax=Int(1e7),gamma=gamma, h_2=1, h_1=1, delta=1.0, delta_2=1.0, L=2^8, mu=[0.1 0.1], hmin=0.2)


# # HI ERN LOGO
# dataall = "data/$(Dates.today())_3D_showcase"
# gamma_0=0.005
# gamma=gamma_0*ones(4,4)
# gamma[2,3]=1.0*gamma_0
# gamma[2,4]=(gamma[2,3]+gamma[3,4])*cospi(1/9)
# run_multilayer(Tmax=Int(1e6),gamma=gamma, h_2=2.0, h_1=10.0, delta=1.0, delta_2=1.0, Lx=2^10,Ly=2^10, mu=[0.1 0.1], hmin=0.1, devicenumber=1, uselogo=true)
