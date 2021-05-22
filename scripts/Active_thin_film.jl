using DrWatson
@quickactivate :Swalbe
using CSV
# Take default values
sys = Swalbe.SysConst_1D()
# Kapital gamma, mobility M, diffusion coefficient D and dumping interval
GG = 0.5
Mob = 1.0
Dif = 0.1
time_dump = 100

"""
    run_active_thin_film()

Simulation with `active` colloids given by a field rho.
"""
function run_active_thin_film(
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