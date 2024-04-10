using Swalbe, Plots, DelimitedFiles, FFTW, CUDA


data = "data/08.03.lbm_vs_ftcs_00"

"""
    run_active_thin_film()

Simulation with `active` colloids given by a field rho.
"""
function run_active_thin_film_old(
    sys::Swalbe.SysConst_1D;
    h₀=1,
    ϵ=0.001,
    θ₀=1/9,  
    verbos=true, 
    dump = time_dump,
    rho = zeros(sys.L),
    differentials = zeros(sys.L, 4), 
    fluid=zeros(sys.Tmax÷dump, sys.L),
    rho_evo=zeros(sys.Tmax÷dump, sys.L),
    Gam = GG,
    mobl = Mob, 
    difu = Dif,
    T=Float64
)
    println("Simulating a active thin film")
    fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p, fluc = Swalbe.Sys(sys, true, T)
    Swalbe.randinterface!(height, h₀, ϵ)
    # Swalbe.sinewave1D!(height, h₀, 1, ϵ, 1)
    rho .= 0.1
    Swalbe.equilibrium!(feq, height, vel)
    ftemp .= feq
    println("Starting the lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(height)
            if verbos
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        # Time update of rho
        Swalbe.update_rho!(rho, fluc, height, dgrad, differentials, M=mobl, D=difu)
        # Film pressure with rho field
        Swalbe.filmpressure!(pressure, height, dgrad, rho, sys.γ, θ₀, sys.n, sys.m, sys.hmin, sys.hcrit, Gamma=Gam)
        Swalbe.∇f!(h∇p, pressure, dgrad, height)
        Swalbe.slippage!(slip, height, vel, sys.δ, sys.μ)
        # HereSwalbe.thermal!(fluc, height, sys.kbt, sys.μ, sys.δ)
        # Here we a force that is like pull of an inclined plane
        F .= h∇p .+ slip # .+ fluc 
        Swalbe.equilibrium!(feq, height, vel)
        Swalbe.BGKandStream!(fout, feq, ftemp, -F)
        Swalbe.moments!(height, vel, fout)
        # Make a snaphot of the configuration
        Swalbe.snapshot!(fluid, height, t, dumping = dump)
        Swalbe.snapshot!(rho_evo, rho, t, dumping = dump)
    end

    return fluid, rho_evo# height, vel, rho
    
end
function  stefans_simulation()
    # Take default values
    sys = Swalbe.SysConst_1D()
    # Kapital gamma, mobility M, diffusion coefficient D and dumping interval
    GG = 0.5
    Mob = 1.0
    Dif = 0.1
    time_dump = 100# Take default values
    sys = Swalbe.SysConst_1D()
    # Kapital gamma, mobility M, diffusion coefficient D and dumping interval
    GG = 0.5
    Mob = 1.0
    Dif = 0.1
    time_dump = 100
   # Now run the actual simulation
   fluid, rho_evo = run_active_thin_film(sys, Gam = 1.0, mobl = 0.1, difu = 0.1)
   # Push the results into a dict for writting to CSV
   df_fluid = Dict()
   df_rho = Dict()
   # Fill the dict with the time steps recorded
   for t in 1:sys.Tmax÷sys.tdump
       df_fluid["h_$(t*time_dump)"] = fluid[t,:]
       df_rho["rho_$(t*time_dump)"] = rho_evo[t,:]
   end
   # Write the desired CSV files
   CSV.write("height_run_Gamma_$(GG)_Dif_$(Dif)_Mob_$(Mob).csv", df_fluid)
   CSV.write("rho_run_Gamma_$(GG)_Dif_$(Dif)_Mob_$(Mob).csv", df_rho)
   # Simulation is done and output is save to files height_...csv and rho_...csv
   println("Sim done, results in height_{...}.csv and rho_{...}.csv") 
end


# function not_active(;Tmax = 10000)
#     sys = Swalbe.SysConst_1D{Float64}(Tmax=Tmax)
#     state = Swalbe.Sys(sys)
#     Swalbe.randinterface!(state.height,1, 0.001)
#     Swalbe.equilibrium!(state)
#     Swalbe.time_loop(sys, state, verbose=true)
#     writedlm("$(data)/run_random_h_0=1,eps=0.001,T=$(Tmax).png", state.height)
#     writedlm("$(data)/run_random_p_0=1,eps=0.001,T=$(Tmax).png", state.pressure)
#     plot(state.height, label="h", title = "run randomn h_0=1, eps=0.001, T=$(Tmax)")
#     savefig("$(data)/run_random_h_0=1,eps=0.001,T=$(Tmax).png")
#     plot(state.pressure, label="p", title = "run randomn h_0=1, eps=0.001, T=$(Tmax)")
#     savefig("$(data)/run_random_p_h_0=1,eps=0.001,T=$(Tmax).png")
#     return nothing
# end

function no_rho(iter;
    h_0=1,
    D=0, M = 0, Γ = 0, K = 0,
    Tmax=10000, ϵ=0.001, plotit= false
)
    sys = Swalbe.SysConstActive_1D{Float64}(D=D, M =M, Γ = Γ , K = K, Tmax=Tmax)
    state = Swalbe.Sys(sys)
    state.rho .= 0
    Swalbe.randinterface!(state.height, h_0,ϵ)
    #println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    println("Starting Lattice Boltzmann time loop")
    Swalbe.time_loop(sys, state, verbose = true)
    writedlm("$(data)/no_rho_height_h_0=$(h_0),eps=$(ϵ)_$(iter).csv", state.height)
    writedlm("$(data)/no_rho_p_h_0=$(h_0),eps=$(ϵ)_$(iter).csv", state.pressure)
    #writedlm("$(data)/active_rho_rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),T=$(Tmax),eps=$(ϵ).csv", state.rho)
    if plotit
        plot(state.height , label = "height", title= "h_0=$(h_0),eps=$(ϵ)")
        #plot!(state.rho, label ="rho")
        savefig("$(data)/no_rho_plot_h_0=$(h_0),T=$(Tmax),eps=$(ϵ)_$(iter).png")
        plot(state.pressure , label = "p", title= "h_0=$(h_0),eps=$(ϵ)")
        #plot!(state.rho, label ="rho")
        savefig("$(data)/no_rho_plot_p_h_0=$(h_0),T=$(Tmax),eps=$(ϵ)_$(iter).png")
    end
end


function active(;
    D=4.5e-5, M = 0.1, Γ = 0.00001, K = 0.0,
    Tmax=100000, tdump=Int(floor(Tmax/1)),
    h= 1,  rho=1, ϵ=0.00001
)
    Plots.gr()
    mkdir(data)
    mkdir("$(data)/figs")
    sys = Swalbe.SysConstActive_1D{Float64}(L =  2^10, D=D, M =M, Γ = Γ , K = K, Tmax=Tmax, tdump=tdump, γ_0 = 0.001)
    state = Swalbe.Sys(sys)
    if rho==0
        state.rho .= 0
    else
        Swalbe.randinterface!(state.rho, rho, ϵ)
    end
    Swalbe.randinterface!(state.height, h,0.0001)
    mass = 0.0
    mass = sum(state.height)
    coloids = 0.0
    coloids  = sum(state.rho)
    println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    println("Starting Lattice Boltzmann time loop for active thin film")
    i=0
    t=0
    # Plots.plot(state.height , label = "height")
    # Plots.plot!(state.rho, label ="rho")
    # Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).png")
    massar = []
    colloidsar = []
    diffar=[]
    append!(diffar, maximum(state.height)-minimum(state.height))
    append!(massar, mass)
    append!(colloidsar, coloids)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            i+=1
            mass = 0.0
            mass = sum(state.height)
            coloids = 0.0
            coloids = sum(state.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.height)
            writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_h_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.γ)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.γ)
            # Plots.plot(state.height , label = "height")
            # p1=Plots.plot!(state.rho, label ="rho", title="t=$(t)")
            # p2=Plots.plot(state.grad_rho, label ="grad rho", title= "colloids = $(coloids)" )
            # p3=Plots.plot(state.grad_h, label ="grad h", title = " mass = $(mass)")
            # p4=Plots.plot(state.γ, label ="gamma")
            # p4=Plots.plot!(state.∇γ, label ="grad gamma")
            # Plots.plot(p1,p2,p3,p4)
            # Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_h=$(h)rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).png")
            append!(massar, mass)
            append!(colloidsar, coloids)
            append!(diffar, maximum(state.height)-minimum(state.height))
        end
        Swalbe.update_rho_LBM!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        Swalbe.∇γ!(state)
        #the sign -3/2 is already in ∇γ! that may be confussing but I want it to be consistent with Stefans code
        state.F .= -state.h∇p .- state.slip .- state.∇γ
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    writedlm("$(data)/active_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", massar)
    writedlm("$(data)/active_colloids_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", colloidsar)
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffar)
    # Plots.plot(massar, label="mass")
    # Plots.plot!(colloidsar, label="colloids", xlabel="$(sys.tdump) t")
    # Plots.savefig("$(data)/figs/active_mass_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
    # Plots.plot(diffar, label="active")
    # Plots.savefig("$(data)/figs/active_diff_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
end



function active_LBM_vs_FTCS(;
    D=4.5e-5, M = 0.1, Γ = 0.00001, K = 0.0,
    Tmax=40000000, tdump=Int(floor(Tmax/100)),
    h= 1,  rho=1, ϵ=0.00001
)
    Plots.gr()
    mkdir(data)
    mkdir("$(data)/figs")
    sys = Swalbe.SysConstActive_1D{Float64}(L =  2^10, D=D, M =M, Γ = Γ , K = K, Tmax=Tmax, tdump=tdump, γ_0 = 0.001)
    state = Swalbe.Sys(sys)
    stateftcs = Swalbe.Sys(sys)
    if rho==0
        state.rho .= 0
    else
        Swalbe.randinterface!(state.rho, rho, ϵ)
    end
    Swalbe.randinterface!(state.height, h,0.0001)
    mass = 0.0
    mass = sum(state.height)
    coloids = 0.0
    coloids  = sum(state.rho)
    stateftcs.rho.=state.rho
    stateftcs.height.=state.height
    massftcs = 0.0
    massftcs = sum(stateftcs.height)
    coloidsftcs = 0.0
    coloidsftcs  = sum(stateftcs.rho)
    println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    Swalbe.equilibrium!(stateftcs)
    println("Starting Lattice Boltzmann time loop for active thin film")
    i=0
    t=0
    # Plots.plot(state.height , label = "height")
    # Plots.plot!(state.rho, label ="rho")
    # Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).png")
    massar = []
    colloidsar = []
    diffar=[]
    massarftcs = []
    colloidsarftcs = []
    diffarftcs=[]
    append!(diffar, maximum(state.height)-minimum(state.height))
    append!(massar, mass)
    append!(colloidsar, coloids)
    append!(diffarftcs, maximum(stateftcs.height)-minimum(stateftcs.height))
    append!(massarftcs, massftcs)
    append!(colloidsarftcs, coloidsftcs)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            i+=1
            mass = 0.0
            mass = sum(state.height)
            coloids = 0.0
            coloids = sum(state.rho)
            massftcs = 0.0
            massftcs = sum(stateftcs.height)
            coloidsftcs = 0.0
            coloidsftcs  = sum(stateftcs.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.height)
            writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_h_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.γ)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.γ)
            writedlm("$(data)/$(to_4_digits(i))_ftcs_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateftcs.height)
            writedlm("$(data)/$(to_4_digits(i))_ftcs_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateftcs.rho)
            writedlm("$(data)/$(to_4_digits(i))_ftcs_active_grad_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateftcs.rho)
            writedlm("$(data)/$(to_4_digits(i))_ftcs_active_grad_h_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateftcs.rho)
            writedlm("$(data)/$(to_4_digits(i))_ftcs_active_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateftcs.γ)
            writedlm("$(data)/$(to_4_digits(i))_ftcs_active_grad_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateftcs.γ)
            Plots.plot(state.height , label = "height lbm")
            p1=Plots.plot!(state.rho, label ="rho lbm", title="t=$(t)")
            p1=Plots.plot!(stateftcs.height , label = "height ftcs")
            p1=Plots.plot!(stateftcs.rho, label ="rho ftcs", title="t=$(t)")
            p2=Plots.plot(state.grad_rho, label ="grad rho lbm", title= "colloids = $(coloids)" )
            p2=Plots.plot!(stateftcs.grad_rho, label ="grad rho ftcs", title= "colloids ftcs= $(coloidsftcs)" )
            p3=Plots.plot(state.grad_h, label ="grad h lbm", title = " mass = $(mass)")
            p3=Plots.plot!(stateftcs.grad_h, label ="grad h ftcs", title = " mass ftcs= $(massftcs)")
            p4=Plots.plot(state.γ, label ="gamma lbm")
            p4=Plots.plot!(state.∇γ, label ="grad gamma lbm")
            p4=Plots.plot!(stateftcs.γ, label ="gamma ftcs")
            p4=Plots.plot!(stateftcs.∇γ, label ="grad gamma ftcs")
            Plots.plot(p1,p2,p3,p4)
            Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_h=$(h)rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).png")
            append!(massar, mass)
            append!(colloidsar, coloids)
            append!(diffar, maximum(state.height)-minimum(state.height))
            append!(diffarftcs, maximum(stateftcs.height)-minimum(stateftcs.height))
            append!(massarftcs, massftcs)
            append!(colloidsarftcs, coloidsftcs)
        end
        Swalbe.update_rho_LBM!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        Swalbe.∇γ!(state)
        #the sign -3/2 is already in ∇γ! that may be confussing but I want it to be consistent with Stefans code
        state.F .= -state.h∇p .- state.slip .- state.∇γ
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        Swalbe.update_rho!(stateftcs, sys)
        Swalbe.filmpressure!(stateftcs, sys)
        Swalbe.h∇p!(stateftcs)
        Swalbe.slippage!(stateftcs, sys)
        Swalbe.∇γ!(stateftcs)
        #the sign -3/2 is already in ∇γ! that may be confussing but I want it to be consistent with Stefans code
        stateftcs.F .= -stateftcs.h∇p .- stateftcs.slip .- stateftcs.∇γ
        Swalbe.equilibrium!(stateftcs)
        Swalbe.BGKandStream!(stateftcs, sys)
        Swalbe.moments!(stateftcs)
    end
    writedlm("$(data)/active_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", massar)
    writedlm("$(data)/active_colloids_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", colloidsar)
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffar)
    writedlm("$(data)/ftcs_active_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", massarftcs)
    writedlm("$(data)/ftcs_active_colloids_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", colloidsarftcs)
    writedlm("$(data)/ftcs_diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffarftcs)
    Plots.plot(massar, label="mass lbm")
    Plots.plot!(colloidsar, label="colloids lbm", xlabel="t/$(sys.tdump)", ylabel="mass")
    Plots.plot!(massarftcs, label="mass ftcs")
    Plots.plot!(colloidsar, label="colloids ftcs", xlabel="t/$(sys.tdump)")
    Plots.savefig("$(data)/figs/active_mass_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
    Plots.plot(diffar, label="lbm")
    Plots.plot!(diffarftcs, label="ftcs", xlabel="t/$(sys.tdump)", ylabel="h_max-h_min")
    Plots.savefig("$(data)/figs/active_diff_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
end

function activeCUDA(;
    D=4.5e-5, M = 0.1, Γ = 0.00001, K = 0.0,
    Tmax=100000, tdump=Int(floor(Tmax/1)),
    h= 1,  rho=1, ϵ=0.00001
)
    println("CUDA")
    # mkdir(data)
    CUDA.allowscalar(false)
    sys = Swalbe.SysConstActive_1D{Float64}(L =  2^10, D=D, M =M, Γ = Γ , K = K, Tmax=Tmax, tdump=tdump, γ_0 = 0.001)
    state = Swalbe.Sys(sys, device="GPU")
    if rho==0
        state.rho .= 0
    else
        Swalbe.Curandinterface!(state.rho, rho, ϵ)
    end
    Swalbe.Curandinterface!(state.height, h,0.0001)
    mass = 0.0
    mass = sum(state.height)
    coloids = 0.0
    coloids  = sum(state.rho)
    println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    println("Starting Lattice Boltzmann time loop for active thin film")
    i=0
    t=0
    massar = []
    colloidsar = []
    diffar=[]
    append!(diffar, maximum(state.height)-minimum(state.height))
    append!(massar, mass)
    append!(colloidsar, coloids)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            i+=1
            mass = 0.0
            mass = sum(state.height)
            coloids = 0.0
            coloids = sum(state.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            # writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.height)
            # writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            # writedlm("$(data)/$(to_4_digits(i))_active_grad_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            # writedlm("$(data)/$(to_4_digits(i))_active_grad_h_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            # writedlm("$(data)/$(to_4_digits(i))_active_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.γ)
            # writedlm("$(data)/$(to_4_digits(i))_active_grad_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.γ)
            append!(massar, mass)
            append!(colloidsar, coloids)
            append!(diffar, maximum(state.height)-minimum(state.height))
        end
        Swalbe.update_rho_LBM!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        Swalbe.∇γ!(state)
        #the sign -3/2 is already in ∇γ! that may be confussing but I want it to be consistent with Stefans code
        state.F .= -state.h∇p .- state.slip .- state.∇γ
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    writedlm("$(data)/active_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", massar)
    writedlm("$(data)/active_colloids_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", colloidsar)
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffar)
end

function minimalistic_active(;
    D=4.5e-7, M = 0.1,  K = 0.0,
    Tmax=50000000, tdump=Int(floor(Tmax/200)),
    h= 1,  rho=1, ϵ=0.00001, γ_0 = 0.01, Γ = 0.01*γ_0*rho, n=3,m=2
)
    mkdir(data)
    mkdir("$(data)/figs")
    sys = Swalbe.SysConstActive_1D{Float64}(L =  2^10, D=D, M =M, Γ = Γ , K = K, Tmax=Tmax, tdump=tdump, γ_0 = γ_0, n=n,m=m)
    state = Swalbe.Sys(sys)
    if rho==0
        state.rho .= 0
    else
        Swalbe.randinterface!(state.rho, rho, ϵ)
    end
    Swalbe.randinterface!(state.height, h,0.0001)
    mass = 0.0
    mass = sum(state.height)
    coloids = 0.0
    coloids  = sum(state.rho)
    println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    println("Starting Lattice Boltzmann time loop for active thin film")
    i=0
    t=0
    massar = []
    colloidsar = []
    diffar=[]
    append!(diffar, maximum(state.height)-minimum(state.height))
    append!(massar, mass)
    append!(colloidsar, coloids)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            i+=1
            mass = 0.0
            mass = sum(state.height)
            coloids = sum(state.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", state.height)
            writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_h_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", state.γ)
            writedlm("$(data)/$(to_4_digits(i))_active_grad_gamma_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", state.γ)
            append!(massar, mass)
            append!(colloidsar, coloids)
            append!(diffar, maximum(state.height)-minimum(state.height))
        end
        if isnan(coloids)
            break
        end
        Swalbe.update_rho_LBM!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        Swalbe.∇γ!(state)
        #the sign -3/2 is already in ∇γ! that may be confussing but I want it to be consistent with Stefans code
        state.F .= -state.h∇p .- state.slip .- state.∇γ
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    writedlm("$(data)/active_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", massar)
    writedlm("$(data)/active_colloids_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv", colloidsar)
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).csv",diffar)
    Plots.plot(massar, label="mass")
    Plots.plot!(colloidsar, label="colloids", xlabel="$(sys.tdump) t")
    Plots.savefig("$(data)/figs/active_mass_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).png")
    Plots.plot(diffar, label="active")
    Plots.savefig("$(data)/figs/active_diff_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ),gamma_0=$(γ_0),n=$(n),m=$(m).png")
end


function compare_active(;
    D=4.5e-5, M = 0.1, Γ = 0.00001, K = 0.0,
    Tmax=10000000, tdump=Int(floor(Tmax/100)),
    h= 1,  rho=1, ϵ=0.00001
)
    mkdir(data)
    mkdir("$(data)/figs")
    Plots.gr()
    sys = Swalbe.SysConstActive_1D{Float64}(L =  h*2^10, D=D, M =M, Γ = Γ ,γ_0 = 0.001, K = K, Tmax=Tmax, tdump=tdump)
    state = Swalbe.Sys(sys)
    sysanti = Swalbe.SysConstActive_1D{Float64}(L =  h*2^10, D=D, M =M, Γ = -Γ ,γ_0 = 0.001, K = K, Tmax=Tmax, tdump=tdump)
    stateanti = Swalbe.Sys(sysanti)
    if rho==0
        state.rho .= 0
    else
        Swalbe.randinterface!(state.rho, rho, ϵ)
    end
    Swalbe.randinterface!(state.height, h,0.0001)
    stateanti.height .= state.height
    stateanti.rho .= state.rho
    mass = 0.0
    mass = sum(state.height)
    coloids = 0.0
    coloids  = sum(state.rho)
    colloidsanti = 0.0
    colloidsanti = sum(stateanti.rho)
    sysn = Swalbe.SysConst_1D{Float64}(L =  h*2^10, Tmax=Tmax, tdump=tdump)
    staten = Swalbe.Sys(sysn)
    staten.height .= state.height
    mass = 0.0
    mass = sum(staten.height)
    #println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    Swalbe.equilibrium!(staten)
    println("Starting Lattice Boltzmann time loop for active thin film")
    i=0
    t=0
    Plots.plot(state.height , label = "height")
    p1=Plots.plot!(state.rho, label ="rho", title="active")
    p2 = Plots.plot(staten.height , label = "height", title="not active")
    Plots.plot(p1,p2)
    Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).png")
    diffar = []
    diffarn = []
    colloidsar = []
    diffaranti = []
    colloidsaranti = []
    append!(diffar, maximum(state.height)-minimum(state.height))
    append!(diffarn, maximum(staten.height)-minimum(staten.height))
    append!(colloidsar, coloids)
    append!(diffaranti, maximum(stateanti.height)-minimum(stateanti.height))
    append!(colloidsaranti, colloidsanti)
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            i+=1
            mass = 0.0
            mass = sum(state.height)
            coloids = 0.0
            coloids = sum(state.rho)
            colloidsanti = 0.0
            colloidsanti = sum(stateanti.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.height)
            writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(-Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateanti.height)
            writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(-Γ),K=$(K),t=$(t),eps=$(ϵ).csv", stateanti.rho)
            writedlm("$(data)/$(to_4_digits(i))_not_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", staten.height)
            # ha=readdlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv")
            # rhoa=readdlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv")
            # hn=readdlm("$(data)/$(to_4_digits(i))_not_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv")
            p1= Plots.plot(state.height , label = "surfactant h")
            p1= Plots.plot!(staten.height , label = "not active h")
            # p1 = Plots.plot!(state.rho, label="surfactant rho", title ="t=$(t)")
            p1= Plots.plot!(stateanti.height , label = "antisurfactant h")
            # p1 = Plots.plot!(stateanti.rho, label="antisurfactant rho", title ="t=$(t)")
            # p1= Plots.plot(ha , label = "active h")
            # p1= Plots.plot!(rhoa , label = "rho")
            # p1 = Plots.plot!(hn, label="not active", title ="t=$(t)")
            fourier = rfft(state.height)
            fouriern = rfft(staten.height)
            fourieranti = rfft(stateanti.height)
            # fourier = rfft(ha)
            # fouriern = rfft(hn)
            @. fourier = ifelse(abs(fourier)==0, NaN, fourier)
            @. fourieranti = ifelse(abs(fourieranti)==0, NaN, fourieranti)
            @. fouriern = ifelse(abs(fouriern)==0, NaN, fouriern)
            writedlm("$(data)/$(to_4_digits(i))_active_fourier_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", fourier)
            writedlm("$(data)/$(to_4_digits(i))_active_fourier_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(-Γ),K=$(K),t=$(t),eps=$(ϵ).csv", fourieranti)
            writedlm("$(data)/$(to_4_digits(i))_not_active_fourier_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).csv", fouriern)
            p2 = Plots.plot(abs.(fourier), xlabel="q", ylabel="σ",label = "surfactant", title="FFT",scale=:log10, legend=:bottomleft)
            p2 = Plots.plot!(abs.(fouriern), xlabel="q", ylabel="σ",label = "not active", title="FFT",scale=:log10, legend=:bottomleft)
            p2 = Plots.plot!(abs.(fourieranti), xlabel="q", ylabel="σ",label = "antisurfactant", title="FFT",scale=:log10, legend=:bottomleft)
            Plots.plot(p1,p2)
            Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),eps=$(ϵ).png")
            append!(diffar, maximum(state.height)-minimum(state.height))
            append!(diffarn, maximum(staten.height)-minimum(staten.height))
            append!(colloidsar, coloids)
            append!(diffaranti, maximum(stateanti.height)-minimum(stateanti.height))
            append!(colloidsaranti, colloidsanti)
        end
        Swalbe.update_rho_LBM!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        Swalbe.∇γ!(state)
        #the signs are ∇γ gives minus then in the following line we make it a plus and the you have ∂_t h+ div (hu)=0 bringing div(hu) to the left side will finally give a minus as suposed to. 
        state.F .= -state.h∇p .- state.slip .- state.∇γ 
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        Swalbe.update_rho_LBM!(stateanti, sysanti)
        Swalbe.filmpressure!(stateanti, sysanti)
        Swalbe.h∇p!(stateanti)
        Swalbe.slippage!(stateanti, sysanti)
        Swalbe.∇γ!(stateanti)
        #the signs are ∇γ gives minus then in the following line we make it a plus and the you have ∂_t h+ div (hu)=0 bringing div(hu) to the left side will finally give a minus as suposed to. 
        stateanti.F .= -stateanti.h∇p .- stateanti.slip .- stateanti.∇γ 
        Swalbe.equilibrium!(stateanti)
        Swalbe.BGKandStream!(stateanti, sysanti)
        Swalbe.moments!(stateanti)
        Swalbe.filmpressure!(staten, sysn)
        Swalbe.h∇p!(staten)
        Swalbe.slippage!(staten, sysn)
        staten.F .= -staten.h∇p .- staten.slip
        Swalbe.equilibrium!(staten)
        Swalbe.BGKandStream!(staten, sysn)
        Swalbe.moments!(staten)
    end
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffar)
    writedlm("$(data)/diff_not_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffarn)
    writedlm("$(data)/colloids_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", colloidsar)
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(-Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv",diffaranti)
    writedlm("$(data)/colloids_mass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(-Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv", colloidsaranti)
    # dif=readdlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv")
    # difn=readdlm("$(data)/diff_not_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).csv")
    Plots.plot(diffar, label="surfactant")
    Plots.plot!(diffaranti, label="antisurfactant")
    Plots.plot!(diffarn, label="not active", xlabel="$(sys.tdump) t", ylabel="h_max - h_min", legend=:upleft)
    Plots.savefig("$(data)/figs/active_diff_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
    Plots.plot(diffar, label="surfactant", scale=:log10)
    Plots.plot!(diffaranti, label="antisurfactant", scale=:log10)
    Plots.plot!(diffarn, label="not active", xlabel="$(sys.tdump) t", ylabel="h_max - h_min", scale=:log10, legend=:upleft)
    Plots.savefig("$(data)/figs/active_diff_log_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
    Plots.plot(colloidsar, label="colloids", ylabel="mass", xlabel="$(sys.tdump) t")
    Plots.plot(colloidsaranti, label="colloids anti", ylabel="mass", xlabel="$(sys.tdump) t")
    Plots.savefig("$(data)/figs/colloidsmass_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
    # Plots.plot(dif, label="active", scale=:log10, legend=:topright)
    # Plots.plot!(difn, label="not active", xlabel="$(sys.tdump)  t", ylabel="h_max - h_min", scale=:log10, legend=:topleft)
    # Plots.savefig("$(data)/figs/active_diff_plot_logscales_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),eps=$(ϵ).png")
end


function compare_active_sin(;
    D=0.0001, M = 0.1, Γ = 0.0005, K = 0.0,
    Tmax=40000000, tdump=10^6,
    h= 1.0,  rho=0.05, amp=0.2
)
    Plots.gr()
    sys = Swalbe.SysConstActive_1D{Float64}(L =  h*2^10, D=D, M =M, Γ = Γ , K = K, Tmax=Tmax, tdump=tdump)
    state = Swalbe.Sys(sys)
    if rho==0
        state.rho .= 0
    else
        Swalbe.othersinewave!(state.rho, rho, 5, amp)
    end
    Swalbe.randinterface!(state.height, h,0.0001)
    mass = 0.0
    mass = sum(state.height)
    coloids = 0.0
    coloids  = sum(state.rho)
    sysn = Swalbe.SysConst_1D{Float64}(L =  h*2^10, Tmax=Tmax, tdump=tdump)
    staten = Swalbe.Sys(sysn)
    staten.height .= state.height
    mass = 0.0
    mass = sum(staten.height)
    #println("Initial mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
    Swalbe.equilibrium!(state)
    Swalbe.equilibrium!(staten)
    println("Starting Lattice Boltzmann time loop for active thin film")
    i=0
    t=0
    Plots.plot(state.height , label = "height")
    p1=Plots.plot!(state.rho, label ="rho", title="active")
    p2 = Plots.plot(staten.height , label = "height", title="not active")
    Plots.plot(p1,p2)
    Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).png")
    diffar = []
    diffarn = []
    append!(diffar, maximum(state.height)-minimum(state.height))
    append!(diffarn, maximum(staten.height)-minimum(staten.height))
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            i+=1
            mass = 0.0
            mass = sum(state.height)
            coloids = 0.0
            coloids = sum(state.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv", state.height)
            writedlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv", state.rho)
            writedlm("$(data)/$(to_4_digits(i))_not_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv", staten.height)
            # ha=readdlm("$(data)/$(to_4_digits(i))_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv")
            # rhoa=readdlm("$(data)/$(to_4_digits(i))_active_rho_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv")
            # hn=readdlm("$(data)/$(to_4_digits(i))_not_active_height_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv")
            p1= Plots.plot(state.height , label = "active h")
            p1= Plots.plot!(staten.height , label = "not active h")
            p1 = Plots.plot!(state.rho, label="rho", title ="t=$(t)")
            # p1= Plots.plot(ha , label = "active h")
            # p1= Plots.plot!(rhoa , label = "rho")
            # p1 = Plots.plot!(hn, label="not active", title ="t=$(t)")
            fourier = rfft(state.height)
            fouriern = rfft(staten.height)
            # fourier = rfft(ha)
            # fouriern = rfft(hn)
            @. fourier = ifelse(abs(fourier)==0, NaN, fourier)
            @. fouriern = ifelse(abs(fouriern)==0, NaN, fouriern)
            writedlm("$(data)/$(to_4_digits(i))_active_fourier_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv", fourier)
            writedlm("$(data)/$(to_4_digits(i))_not_active_fourier_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).csv", fouriern)
            p2 = Plots.plot(abs.(fourier), xlabel="q", ylabel="σ",label = "active", title="FFT",scale=:log10, legend=:bottomleft)
            p2 = Plots.plot!(abs.(fouriern), xlabel="q", ylabel="σ",label = "not active", title="FFT",scale=:log10, legend=:bottomleft)
            Plots.plot(p1,p2)
            Plots.savefig("$(data)/figs/$(to_4_digits(i))_active_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp).png")
            append!(diffar, maximum(state.height)-minimum(state.height))
            append!(diffarn, maximum(staten.height)-minimum(staten.height))
        end
        Swalbe.update_rho!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
        Swalbe.filmpressure!(staten, sysn)
        Swalbe.h∇p!(staten)
        Swalbe.slippage!(staten, sysn)
        staten.F .= -staten.h∇p .- staten.slip
        Swalbe.equilibrium!(staten)
        Swalbe.BGKandStream!(staten, sysn)
        Swalbe.moments!(staten)
    end
    writedlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),amp=$(amp).csv",diffar)
    writedlm("$(data)/diff_not_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),amp=$(amp).csv",diffarn)
    # dif=readdlm("$(data)/diff_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),amp=$(amp).csv")
    # difn=readdlm("$(data)/diff_not_active_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),amp=$(amp).csv")
    Plots.plot(diffar, label="active")
    Plots.plot!(diffarn, label="not active", xlabel="$(sys.tdump) t", ylabel="h_max - h_min")
    Plots.savefig("$(data)/figs/active_diff_plot_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),amp=$(amp).png")
    # Plots.plot(dif, label="active", scale=:log10, legend=:topright)
    # Plots.plot!(difn, label="not active", xlabel="$(sys.tdump) t", ylabel="h_max - h_min", scale=:log10, legend=:topleft)
    # Plots.savefig("$(data)/figs/active_diff_plot_logscales_h=$(h),rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(Tmax),amp=$(amp).png")
end


function to_4_digits(i)
    if i <10
        return "000$(i)"
    elseif i<100
        return "00$(i)"
    elseif i<1000
        return "0$(i)"
    else
        return "$(i)"
    end
end



function active_sin(;
    D=0.0001, M = 0.0, Γ = 0.00001, K = 0.0,
    Tmax=10000, amp=0.01, tdump=Int(floor(Tmax/10)),
    h_0 = 1,  rho_0=1, ϵ=0.0
)
    sys = Swalbe.SysConstActive_1D{Float64}(D=D, M =M, Γ = Γ , K = K, Tmax=Tmax, tdump=tdump)
    state = Swalbe.Sys(sys)
    Swalbe.sinewave1D!(state.rho, rho_0, 2, amp,1)
    Swalbe.randinterface!(state.height, h_0,ϵ)
    Swalbe.stability_analysis(sys,state)
    writedlm("$(data)/initial_data_roller_height_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).csv", state.height)
    writedlm("$(data)/initial_data_roller_rho_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).csv", state.rho)
    Plots.plot(state.height, label="height", title = "initial data")
    p1 = Plots.plot!(state.rho, label="rho", title = "initial data")
    #savefig("$(data)/roller_initial_image_h,rho_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).png")
    Swalbe.equilibrium!(state)
    println("Starting Lattice Boltzmann time loop")
    for t in 1:sys.Tmax
        if t % sys.tdump == 0
            mass = 0.0
            mass = sum(state.height)
            coloids = 0.0
            coloids = sum(state.rho)
            println("Time step $t mass is $(round(mass, digits=3)) and coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/roller_height_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp),eps=$(ϵ).csv", state.height)
            writedlm("$(data)/roller_rho_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp),eps=$(ϵ).csv", state.rho)
            writedlm("$(data)/roller_p_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp),eps=$(ϵ).csv", state.pressure)
            p2 = Plots.plot(state.height, label="height", title = "after $(t) time steps")
            p3= Plots.plot(state.rho, label="rho", title = "after $(t) time steps")
            #savefig("$(data)/roller_image_h,rho_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).png")
            p4 = Plots.plot(state.pressure, label="pressure", title = "after $(t) time steps")
            #savefig("$(data)/roller_image_p_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).png")
            Plots.savefig(Plots.plot(p1,p2,p3,p4), "$(data)/roller_image_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),t=$(t),amp=$(amp),eps=$(ϵ).png")
        end
        Swalbe.update_rho!(state, sys)
        Swalbe.filmpressure!(state, sys)
        Swalbe.h∇p!(state)
        Swalbe.slippage!(state, sys)
        state.F .= -state.h∇p .- state.slip
        Swalbe.equilibrium!(state)
        Swalbe.BGKandStream!(state, sys)
        Swalbe.moments!(state)
    end
    writedlm("$(data)/final_roller_height_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).csv", state.height)
    writedlm("$(data)/final_roller_rho_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).csv", state.rho)
    writedlm("$(data)/final_roller_p_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).csv", state.pressure)
    p2 = plot(state.height, label="height", title = "after $(Tmax) time steps")
    p3= plot(state.rho, label="rho", title = "after $(Tmax) time steps")
    #savefig("$(data)/roller_image_h,rho_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).png")
    p4 = plot(state.pressure, label="pressure", title = "after $(Tmax) time steps")
    #savefig("$(data)/roller_image_p_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).png")
    savefig(plot(p1,p2,p3,p4), "$(data)/final_roller_image_h_0=$(h_0),rho_0=$(rho_0),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),amp=$(amp),eps=$(ϵ).png")
end

function active_3D(;
    rho = 0,
    D=0, M = 0, Γ = 0, K = 0,
    Tmax=1000, ϵ=0
)
    sys = Swalbe.SysConstActive{Float64}(D=D, M =M, Γ = Γ , K = K, Tmax=Tmax)
    state = Swalbe.Sys(sys, "CPU")
    if rho==0
        state.rho .= 0
    else
        Swalbe.randinterface!(state.rho, rho, ϵ)
    end
    Swalbe.randinterface!(state.height, 1,0.001)
    Swalbe.equilibrium!(state)
    #println("Starting Lattice Boltzmann time loop")
    Swalbe.time_loop(sys, state, verbose = false)
    if isnan(state.rho[256,256])
        println("NaN rho=$(rho), D=$(D), K=$(K), M=$(M), Γ=$(Γ), eps=$(ϵ)")
    end
    writedlm("$(data)/3D_active_height_rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),epsilon=$(ϵ).csv", state.height)
    writedlm("$(data)/3D_active_rho_rho=$(rho),D=$(D),M=$(M),Gamma=$(Γ),K=$(K),Tmax=$(Tmax),epsilon=$(ϵ).csv", state.rho)
end

function test()
    sys = SysConstActive_1D{Float64}(D=0.5, M = 0.5, Γ = 0.5 , K = 0.5, Tmax=10000)
    result = run_active_thin_film(sys,state->Swalbe.randinterface!(state, 1, 0.001))
    p1 = plot(result[1], color = 1, label = "height") 
    p2 = plot(result[2], color = 2, label = "rho")
    gui(merge_series!(p1,p2))
end


function merge_series!(sp1::Plots.Subplot, sp2::Plots.Subplot)
    append!(sp1.series_list, sp2.series_list)
    Plots.expand_extrema!(sp1[:xaxis], xlims(sp2))
    Plots.expand_extrema!(sp1[:yaxis], ylims(sp2))
    Plots.expand_extrema!(sp1[:zaxis], zlims(sp2))
    return sp1
end

function merge_series!(plt, plts...)
    for (i, sp) in enumerate(plt.subplots)
        for other_plt in plts
            if i in eachindex(other_plt.subplots)
                merge_series!(sp, other_plt[i])
            end
        end
    end
    return plt
end

function all_3D_simualtion()
    for K in [0,0.5], M in [0,0.5], Γ in [0,0.5]
        try
            active_3D(K=K, M=M, Γ=Γ, )
        catch
            println("failed rho=0, D=0, K=$(K), M=$(M), Γ=$(Γ), eps=0")
        end
        try
            active_3D(K=K, M=M, Γ=Γ, rho=0.1)
        catch
            println("failed rho=1, D=0, K=$(K), M=$(M), Γ=$(Γ), eps=0")
        end
        try
            active_3D(K=K, M=M, Γ=Γ, rho=0.1, D=0.5)
        catch
            println("failed rho=1, D=0.5, K=$(K), M=$(M), Γ=$(Γ), eps=0")
        end
        try
            active_3D(K=K, M=M, Γ=Γ, rho=0.1, ϵ=0.001)
        catch
            println("failed rho=1, D=0, K=$(K), M=$(M), Γ=$(Γ), eps=0.01")
        end
        try
            active_3D(K=K, M=M, Γ=Γ, rho=0.1, ϵ=0.001, D=0.5)
        catch
            println("failed rho=1, D=0.5, K=$(K), M=$(M), Γ=$(Γ), eps=0.01")
        end
    end
    return nothing
end

function all_2D_simulations()
    for K in [0,0.01], M in [0,0.01], Γ in [0,0.01]
        try
            active(K=K, M=M, Γ=Γ, )
        catch
            println("failed rho=0, D=0, K=$(K), M=$(M), Γ=$(Γ), eps=0")
        end
        try
            active(K=K, M=M, Γ=Γ, rho=0.1)
        catch
            println("failed rho=0.1, D=0, K=$(K), M=$(M), Γ=$(Γ), eps=0")
        end
        try
            active(K=K, M=M, Γ=Γ, rho=0.1, D=0.5)
        catch
            println("failed rho=0.1, D=0.5, K=$(K), M=$(M), Γ=$(Γ), eps=0")
        end
        try
            active(K=K, M=M, Γ=Γ, rho=0.1, ϵ=0.001)
        catch
            println("failed rho=0.1, D=0, K=$(K), M=$(M), Γ=$(Γ), eps=0.01")
        end
        try
            active(K=K, M=M, Γ=Γ, rho=0.1, ϵ=0.001, D=0.5)
        catch
            println("failed rho=0.1, D=0.5, K=$(K), M=$(M), Γ=$(Γ), eps=0.01")
        end
    end
    return nothing
end

function test_update_rho(; h_0 = 1, rho_0=1, Tmax = 10000000, L = 256, D=0.25, amp = 0.5,update_h=false, tdump = Int(floor(Tmax/10)))
    sys = Swalbe.SysConstActive_1D(D=D, L=L)
    state = Swalbe.Sys(sys)
    Swalbe.othersinewave!(state.height , 1, 2, amp)
    # state.height .= h_0
    state.rho .= rho_0
    # Swalbe.sinewave1D!(state.height , 1, 2, amp,1)
    #Swalbe.othersinewave!(state.rho , 1, 1, amp)
    writedlm("$(data)/test_update_rho$(update_h ? "_free_h" : "")_height_h_0=$(h_0),rho_0=$(rho_0),t=0_D=$(D)_L=$(L)_amp=$(amp).csv", state.height)
    writedlm("$(data)/test_update_rho$(update_h ? "_free_h" : "")_rho_h_0=$(h_0),rho_0=$(rho_0),t=$(0)_D$(D)_L=$(L)_$(amp).csv", state.rho)
    plot(state.height, label="height", title = "t=0")
    plot!(state.rho, label="rho")
    i=0
    Swalbe.rho_equilibrium!(state)
    savefig("$(data)/figs/$(to_4_digits(i))_test_update_rho$(update_h ? "_free_h" : "")_h_0=$(h_0),rho_0=$(rho_0),t=0_D=$(D)_L=$(L)_amp=$(amp).png")
    for t in 1:Tmax
        Swalbe.update_rho!(state,sys)
        if t%tdump==0
            i+=1
            coloids = 0.0
            coloids  = sum(state.rho)
            mass = 0.0
            mass = sum(state.height)
            println("Time step is $(t) mass is $(mass) coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/test_update_rho$(update_h ? "_free_h" : "")_rho_h_0=$(h_0),rho_0=$(rho_0),t=$(t)_D$(D)_L=$(L)_$(amp).csv", state.rho)
            Plots.plot(state.rho, label="rho", title="t=$(t)")
            Plots.plot!(state.grad_rho, label="grad rho", title="t=$(t)")
            Plots.plot!(state.lap_rho, label="lap rho", title="t=$(t)")
            Plots.plot!(state.grad_h, label="grad_h", title="t=$(t)")
            Plots.plot!(state.lap_h, label="lap_h", title="t=$(t)")
            if update_h
                writedlm("$(data)/test_update_rho$(update_h ? "_free_h" : "")_height_rho_h_0=$(h_0),rho_0=$(rho_0),t=$(t)_D$(D)_L=$(L)_$(amp).csv", state.height)
                Plots.plot!(state.height, label="h", title="t=$(t)")
            end
            Plots.savefig("$(data)/figs/$(to_4_digits(i))test_update_rho$(update_h ? "_free_h" : "")_h_0=$(h_0),rho_0=$(rho_0),t=$(t)_D=$(D)_L=$(L)_amp=$(amp).png")
        end
        if update_h
            Swalbe.filmpressure!(state, sys)
            Swalbe.h∇p!(state)
            Swalbe.slippage!(state, sys)
            state.F .= -state.h∇p .- state.slip
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state, sys)
            Swalbe.moments!(state)
        end
    end
end

function test_update_rho_drop(; h_0 = 1, rho_0=1, Tmax = 10000, L = 248, D=1e-5, amp = 0.1)
    sys = Swalbe.SysConstActive_1D(D=D, L=L)
    state = Swalbe.Sys(sys)
    Swalbe.singledroplet!(state , 40, 1/9, 124)
    state.rho .= rho_0
    tdump= Int(floor(Tmax/10))
    for t in 1:Tmax
        Swalbe.update_rho!(state,sys)
        if t%tdump==0
            coloids = 0.0
            coloids  = sum(state.rho)
            println("Time step is $(t) coloid-mass is $(round(coloids, digits=3))")
        end
    end
    writedlm("$(data)/test_update_rho_drop_height_h_0=$(h_0),rho_0=$(rho_0),Tmax=$(Tmax)_D=$(D)_L=$(L)_amp=$(amp).csv", state.height)
    writedlm("$(data)/test_update_rho_drop_rho_h_0=$(h_0),rho_0=$(rho_0)Tmax=$(Tmax)_D$(D)_L=$(L)_$(amp).csv", state.rho)
    plot(state.height, label="height", title = "h_0=$(h_0),rho_0=$(rho_0)Tmax=$(Tmax)_D=$(D)")
    plot!(state.rho, label="rho")
    savefig("$(data)/figs/THIS_h_0=$(h_0),rho_0=$(rho_0),Tmax=$(Tmax)_D=$(D)_L=$(L)_$(amp).png")
end




function create_droplet(;L = 248, r=40, theta = 1/9, center=124, Tmax=1000000)
    sys = Swalbe.SysConst_1D(L=L, Tmax=Tmax, θ=theta)
    state = Swalbe.Sys(sys)
    Swalbe.singledroplet!(state.height , r, theta, center)
    Swalbe.time_loop(sys, state, verbose=true)
    return state.height
end

function test_update_rho_droplet(L = 248, r=40, theta = 1/9, center=124, Tmaxdrop=1000000, Tmax = 1000, D=1e-5, K=0.0, M=0.0, Γ=0.00,rho_0=1,r_rho=30,theta_rho=1/9, update_h=false)
    tdump= 100
    sys = Swalbe.SysConstActive_1D(L=L, Tmax=Tmax, K=K, M=M, Γ=Γ)
    state = Swalbe.Sys(sys)
    state.height .= readdlm("$(data)/relaxed_drop_r=$(r),theta=$(theta),center=$(center),Tmax=$(Tmaxdrop).csv")
    Swalbe.singledroplet!(state.rho , r_rho, theta_rho, center)
    t=0
    writedlm("$(data)/test_relaxed_drop_rho_drop_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),r_rho=$(r_rho),theta_rho=$(theta_rho).csv", state.rho)
    plot(state.height, label="height", title = "t=$(t)")
    plot!(state.rho, label="rho")
    savefig("$(data)/test_relaxed_drop_rho_drop_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),r_rho=$(r_rho),theta_rho=$(theta_rho).png")
    for t in 1:Tmax
        Swalbe.update_rho!(state,sys)
        if t%tdump==0
            coloids = 0.0
            coloids  = sum(state.rho)
            println("Time step is $(t) coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/test_relaxed_drop_rho_drop_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),r_rho=$(r_rho),theta_rho=$(theta_rho).csv", state.rho)
            plot(state.height, label="height", title = "t=$(t)")
            plot!(state.rho, label="rho")
            savefig("$(data)/test_relaxed_drop_rho_drop_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),r_rho=$(r_rho),theta_rho=$(theta_rho).png")
        end
        if update_h
            Swalbe.filmpressure!(state, sys)
            Swalbe.h∇p!(state)
            Swalbe.slippage!(state, sys)
            state.F .= -state.h∇p .- state.slip
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state, sys)
            Swalbe.moments!(state)
        end
    end
end


function test_update_rho_double_droplet(;L = 248, r=40, theta = 1/9, center=124, Tmaxdrop=1000000, Tmax = 1000, D=1e-5, K=0.0, M=0.0, Γ=0.00,alpha = 1/2, update_h=false)
    sys = Swalbe.SysConstActive_1D(L=L, Tmax=Tmax, K=K, M=M, Γ=Γ)
    state = Swalbe.Sys(sys)
    height = create_droplet(L = L, r=r, theta = theta, center=center, Tmax=Tmaxdrop)
    #rho= create_droplet(L = L, r=r_rho, theta = theta_rho, center=center, Tmax=Tmaxdrop)
    state.height .= height
    state.rho .= alpha .* height
    t=0
    writedlm("$(data)/test_relaxed_double_drop$(update_h ? "_free_h" : "")_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),alpha=$(alpha).csv", state.rho)
    plot(state.height, label="height", title = "t=$(t)")
    plot!(state.rho, label="rho")
    savefig("$(data)/test_relaxed_drop_double_drop$(update_h ? "_free_h" : "")_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),alpha=$(alpha).png")
    for t in 1:Tmax
        Swalbe.update_rho!(state,sys)
        if t%sys.tdump==0
            coloids = 0.0
            coloids  = sum(state.rho)
            println("Time step is $(t) coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/test_relaxed_drop_double_drop$(update_h ? "_free_h" : "")_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),alpha=$(alpha).csv", state.rho)
            plot(state.height, label="height", title = "t=$(t)")
            plot!(state.rho, label="rho")
            savefig("$(data)/test_relaxed_drop_double_drop$(update_h ? "_free_h" : "")_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),alpha=$(alpha).png")
        end
        if update_h
            Swalbe.filmpressure!(state, sys)
            Swalbe.h∇p!(state)
            Swalbe.slippage!(state, sys)
            state.F .= -state.h∇p .- state.slip
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state, sys)
            Swalbe.moments!(state)
        end
    end
end


function test_update_rho_droplet_line(L = 248, r=40, theta = 1/9, center=124, Tmaxdrop=1000000, Tmax = 1000, D=1e-5, K=0.0, M=0.0, Γ=0.00,rho_0=1,eps=0.001, update_h=false)
    tdump=Tmax/10
    sys = Swalbe.SysConstActive_1D(L=L, Tmax=Tmax, K=K, M=M, Γ=Γ)
    state = Swalbe.Sys(sys)
    state.height .= readdlm("$(data)/relaxed_drop_r=$(r),theta=$(theta),center=$(center),Tmax=$(Tmaxdrop).csv")
    #Swalbe.singledroplet!(state.rho , r_rho, theta_rho, center)
    Swalbe.randinterface!(state.rho, 1, eps)
    t=0
    writedlm("$(data)/test_relaxed_drop_line_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),eps=$(eps).csv", state.rho)
    plot(state.height, label="height", title = "t=$(t)")
    plot!(state.rho, label="rho")
    savefig("$(data)/test_relaxed_drop_line_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),eps=$(eps).png")
    for t in 1:Tmax
        Swalbe.update_rho!(state,sys)
        if t%tdump==0
            coloids = 0.0
            coloids  = sum(state.rho)
            println("Time step is $(t) coloid-mass is $(round(coloids, digits=3))")
            writedlm("$(data)/test_relaxed_drop_line_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),,eps=$(eps).csv", state.rho)
            plot(state.height, label="height", title = "t=$(t)")
            plot!(state.rho, label="rho")
            savefig("$(data)/test_relaxed_drop_line_r=$(r),theta=$(theta),center=$(center),t=$(t),D=$(D),K=$(K),M=$(M),Γ=$(Γ),rho_0=$(rho_0),,eps=$(eps).png")
        end
        if update_h
            Swalbe.filmpressure!(state, sys)
            Swalbe.h∇p!(state)
            Swalbe.slippage!(state, sys)
            state.F .= -state.h∇p .- state.slip
            Swalbe.equilibrium!(state)
            Swalbe.BGKandStream!(state, sys)
            Swalbe.moments!(state)
        end
    end
end




# active(rho = 1,
# D=0, M = 0, Γ = 0.0001, K = 0,
# Tmax=100000, ϵ=0.0001, plotit = true)

t_0 = 8.6* 1e6
Tmax = Int(ceil(3* (t_0)))

#create_droplet()

#test_update_rho_droplet_line()

#test_update_rho(Tmax=10000, tdump = 50)

#test_update_rho_double_droplet(Tmax = Tmax, update_h=true)


# for h in 1:10
#     active(h=h,Tmax = t_0*3*h , tdump=1e6*h, rho=0.2*h, D=1e-3)
# end

# for D in [5*1e-3, 1e-3, 5*1e-4, 1e-4, 5*1e-5, 1e-5, 5*1e-6, 1e-6 ]
#     compare_active(rho=1, D=D)
# end

# compare_active(rho=1, D=0.001,Tmax =10000000, tdump = 1000000)

# for Γ in [0, 0.001, 0,0005, 0,0001]
#     for M in [0, 1, 0.1,0.01,0.001, 0.0001]
#         for sM in [-1,1]
#             for sG in [-1,1]
#                 active(Γ=sG*Γ, M = sM* M )
#             end
#         end 
#     end
# end 


# active(Γ=0.0001, M=1, D=5e-5)

#compare_active()

#active()

# for fac in  [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9]
#     minimalistic_active(D=4.5*fac)
# end
# minimalistic_active(Γ=0)


# minimalistic_active(D=1e-5)
# minimalistic_active(D=1e-7)
# minimalistic_active(D=1e-9)
# minimalistic_active(D=1e-2)
# minimalistic_active(Γ=0)

active_LBM_vs_FTCS()

#active(h=3,Tmax = 5*Tmax, tdump=Tmax/2, rho=0.05, D=1e-3)




#no_rho("02", h_0=1, ϵ= 0.001, Tmax = Tmax,  plotit = true)

#not_active(Tmax=Tmax)


#test_update_rho()


            