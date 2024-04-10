#Parent type for all LBM states 
    abstract type LBM_state end
    # Derived type for 2D LBM states
    abstract type LBM_state_2D <: LBM_state end
    
    abstract type Expanded_2D <: LBM_state_2D end
    # Derived type for 1D LBM states
    abstract type LBM_state_1D <: LBM_state end
    
    abstract type Expanded_1D <: LBM_state_1D end

    abstract type Active_1D <: LBM_state_1D end
    
    # Parent type for system specific constants
    abstract type Consts end
    

"""
    SysConst{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.

# Arguments

- `Lx :: Int`: Length of the lattice sides in x-direction
- `Ly :: Int`: Length of the lattice sides in y-direction
- `Tmax :: Int`: Number of lattice Boltzmann time iterations
- `tdump :: Int`: Dumping interval for e.g. data output
- `τ :: T`: BGK relaxation rate 
- `cₛ :: T`: Lattice speed of sound, every physical velocity needs to be smaller than this! 
- `μ :: T`: Kinematic fluid viscosity
- `δ :: T`: Slip length, defines how far the **no-slip** condition is interpolated into the substrate
- `kbt :: T`: Thermal energy of the film, works with small values ≈ 10^(-7)
- `γ :: T`: Surface tension
- `n :: Int`: Greater exponent of the two used for the powerlaw in the disjoining pressure
- `m :: Int`: Smaller exponent of the two used for the powerlaw in the disjoining pressure
- `hmin :: T`: Height value at which the disjoining pressure functional vanishes
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term
- `θ :: T`: Contact angle in multiples of π
- `g :: T`: gravitational acceleration, usually neglected in thin film simulations

"""
Base.@kwdef struct SysConst{T}
    # Lattice
    Lx :: Int = 256
    Ly :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    # Collision related
    τ :: T = 1.0
    cₛ :: T = 1/sqrt(3.0)
    μ :: T =  cₛ^2 *( τ - 0.5)
    # Force related
    δ :: T = 1.0
    kbt :: T = 0.0
    γ :: T = 0.01
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    θ :: T = 1/9
    g :: T = 0.0
end



"""
SysConstActive{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on. With extra constants for active matter model

# Arguments

- `Lx :: Int`: Length of the lattice sides in x-direction
- `Ly :: Int`: Length of the lattice sides in y-direction
- `Tmax :: Int`: Number of lattice Boltzmann time iterations
- `tdump :: Int`: Dumping interval for e.g. data output
- `τ :: T`: BGK relaxation rate 
- `cₛ :: T`: Lattice speed of sound, every physical velocity needs to be smaller than this! 
- `μ :: T`: Kinematic fluid viscosity
- `δ :: T`: Slip length, defines how far the **no-slip** condition is interpolated into the substrate
- `kbt :: T`: Thermal energy of the film, works with small values ≈ 10^(-7)
- `γ_0 :: T`: initial Surface tension
- `n :: Int`: Greater exponent of the two used for the powerlaw in the disjoining pressure
- `m :: Int`: Smaller exponent of the two used for the powerlaw in the disjoining pressure
- `hmin :: T`: Height value at which the disjoining pressure functional vanishes
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term
- `θ_0 :: T`: initial Contact angle in multiples of π
- `g :: T`: gravitational acceleration, usually neglected in thin film simulations
- `D :: T = 1.0` : TODO
- `M :: T = 0.0` : TODO
- `Γ :: T = 0.0` : TODO
- `K :: T = 0.0` : TODO
- `activemode :: String` : choose from \"surface_tension\", \"contact_angle\" or "" to either calculate γ from ρ 
then θ by Youngs law or θ from ρ and the γ by Youngs law or both from ρ
"""
Base.@kwdef struct SysConstActive{T}
    # Lattice
    Lx :: Int = 256
    Ly :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    # Collision related
    τ :: T = 1.0
    cₛ :: T = 1/sqrt(3.0)
    μ :: T =  cₛ^2 *( τ - 0.5)
    # Force related
    δ :: T = 1.0
    kbt :: T = 0.0
    γ_0 :: T = 0.01
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    θ_0 :: T = 1/9
    g :: T = 0.0
    #active matter realted
    D :: T = 1.0
    Γ :: T = 0.0
    GammaB :: T = 0.0
    rho_BRes_sigma_B_down :: T = 0.0
end

"""
    SysConst_1D{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.

# Arguments

- `L :: Int`: Number of lattice points
- `Tmax :: Int`: Number of lattice Boltzmann time iterations
- `tdump :: Int`: Dumping interval for e.g. data output
- `τ :: T`: BGK relaxation rate 
- `cₛ :: T`: Lattice speed of sound, every physical velocity needs to be smaller than this! 
- `μ :: T`: Kinematic fluid viscosity
- `δ :: T`: Slip length, defines how far the **no-slip** condition is interpolated into the substrate
- `kbt :: T`: Thermal energy of the film, works with small values ≈ 10^(-7)
- `γ :: T`: Surface tension
- `n :: Int`: Greater exponent of the two used for the powerlaw in the disjoining pressure
- `m :: Int`: Smaller exponent of the two used for the powerlaw in the disjoining pressure
- `hmin :: T`: Height value at which the disjoining pressure functional vanishes
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term
- `θ :: T`: Contact angle in multiples of π
- `g :: T`: gravitational acceleration, usually neglected in thin film simulations
"""
Base.@kwdef struct SysConst_1D{T}
    # Lattice
    L :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    # Collision related
    τ :: T = 1.0
    cₛ :: T = 1/sqrt(3.0)
    μ :: T =  cₛ^2 *( τ - 0.5)
    # Force related
    δ :: T = 1.0
    kbt :: T = 0.0
    γ :: T = 0.01
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    θ :: T = 1/9
    g :: T = 0.0
end


"""
    SysConstActive_1D{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on. With extra constants for active matter model.

# Arguments

- `L :: Int`: Number of lattice points
- `Tmax :: Int`: Number of lattice Boltzmann time iterations
- `tdump :: Int`: Dumping interval for e.g. data output
- `τ :: T`: BGK relaxation rate 
- `cₛ :: T`: Lattice speed of sound, every physical velocity needs to be smaller than this! 
- `μ :: T`: Kinematic fluid viscosity
- `δ :: T`: Slip length, defines how far the **no-slip** condition is interpolated into the substrate
- `kbt :: T`: Thermal energy of the film, works with small values ≈ 10^(-7)
- `γ_0 :: T`: Surface tension
- `n :: Int`: Greater exponent of the two used for the powerlaw in the disjoining pressure
- `m :: Int`: Smaller exponent of the two used for the powerlaw in the disjoining pressure
- `hmin :: T`: Height value at which the disjoining pressure functional vanishes
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term
- `θ_0 :: T`: Contact angle in multiples of π
- `g :: T`: gravitational acceleration, usually neglected in thin film simulations
-  `D :: T`: Diffusion constant of the catalyst
- `R :: T`: catalyst radius, not used
- `Γ :: T`: Surface tension change parameter
- `adv_limit :: T`: no idea anymore, Probably depricated
- `rho_crit :: T`: parameter controling the precursor/disjoining pressure of rho see `rho_pressure!`
- `rho_n :: Int`: power of the precursor/disjoining pressure of rho see `rho_pressure!` deprecated
- `rho_pow :: Int`: also from the `rho_pressure` and deprecated
- `rho_max :: T`: No idea anymore, Probably deprecated
- `alpha :: T`: vertical distribution of catalyst
- `A_2 :: T`: second order diffusion correction, not used very often
- `weight_center :: T`: moving averadge smoothing weight
- `production_rate :: T`: production rate
- `sink_rate_bulk :: T`: deprecated reaction model
- `sink_rate_evap :: T`: deprecated reaction model
- `source_rate_reactant :: T`: deprecated reaction model
- `sigma_A_up :: T`: volatility rate of product
- `sigma_B_up :: T`: volatility rate of reactant
- `rho_BRes_sigma_B_down :: T`: source rate of reactant
- `D_Product :: T`: Diffusion coefficient of product
- `D_B :: T`:  Diffusion coefficient of reactant
- `rho_BRes :: T`: reservoir concentration of reactant
- `GammaB :: T`: surface tension change by reactant
- ` hcrit2 :: T`: extra dijoining pressure parameter for `rho_pressure!`, `rho_A_pressure!`, `rho_B_pressure!`
- `b_B :: T`: controling `rho_B_pressure!`
- `b_A :: T`: controling `rho_A_pressure!`
- prescursor_nabla_F :: T`: value of `nabla_F` in the precursor, helps prevent unwanted diffusion into the precursor film
"""
Base.@kwdef struct SysConstActive_1D{T}
    # Lattice
    L :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    # Collision related
    τ :: T = 1.0
    cₛ :: T = 1/sqrt(3.0)
    μ :: T =  cₛ^2 *( τ - 0.5)
    # Force related
    δ :: T = 1.0
    kbt :: T = 0.0
    γ_0 :: T = 0.01
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    θ_0 :: T = 1/9
    g :: T = 0.0
    #active matter realted
    D :: T = 0.01
    R :: T = 0.0
    Γ :: T = 0.0
    adv_limit :: T = 0.0
    rho_crit :: T = 0.05
    rho_n :: Int = 5
    rho_pow :: Int = 8
    rho_max :: T = 3.0
    alpha :: T = 0.0
    A_2 :: T = 0.0
    weight_center :: T = 1.0
    production_rate :: T =0.1
    sink_rate_bulk :: T=0.1
    sink_rate_evap :: T=0.1
    source_rate_reactant :: T=0.1
    sigma_A_up :: T =0.1
    sigma_B_up :: T =0.1
    rho_BRes_sigma_B_down :: T =0.1
    D_Product :: T=6e-4
    D_B :: T = D_Product
    #to be adopted for the four field method for both rho_A and rho_B
   rho_BRes :: T = 0.0
    GammaB :: T = 0.0
    hcrit2 :: T =-0.04
    b_B :: T = (hmin-hcrit+hcrit2 - rho_BRes_sigma_B_down/((hmin - hcrit + hcrit2-rho_crit)*production_rate+sigma_B_up/(hmin-hcrit))) <= 0.0 || isnan(hmin-hcrit+hcrit2 - rho_BRes_sigma_B_down/((hmin-hcrit+hcrit2-rho_crit)*production_rate+sigma_B_up/(hmin-hcrit))) ? 0.01 : (hmin-hcrit+hcrit2 - rho_BRes_sigma_B_down/((hmin-hcrit+hcrit2-rho_crit)*production_rate+sigma_B_up/(hmin-hcrit))) 
    b_A :: T = (hmin-hcrit+hcrit2- production_rate*(hmin-hcrit+hcrit2-rho_crit)*(hmin-b_B)*(hmin-hcrit)^2/sigma_A_up) <=0.0 || isnan(hmin-hcrit+hcrit2- production_rate*(hmin-hcrit+hcrit2-rho_crit)*(hmin-hcrit+hcrit2-b_B)*(hmin-hcrit)^2/sigma_A_up) ? 0.01 : (hmin-hcrit+hcrit2- production_rate*(hmin-hcrit+hcrit2-rho_crit)*(hmin-hcrit+hcrit2-b_B)*(hmin-hcrit)^2/sigma_A_up)
    prescursor_nabla_F :: T = 1/(hmin - hcrit)
end

"""
    State{T, N}

Data structure that stores all arrays for a given simulation.

# Arguments

- `fout :: Array{T,N}`: Output distribution function
- `ftemp :: Array{T,N}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: Array{T,N}`: Equilibrium distribution function
- `height :: Matrix{T}`: Field that stores the scalar height values
- `velx :: Matrix{T}`: Field that stores the x-component of the velocity vector
- `vely :: Matrix{T}`: Field that stores the y-component of the velocity vector
- `vsq :: Matrix{T}`: Field that stores the velocity squared, used in `equilibrium!`
- `pressure :: Matrix{T}`: Pressure distribution computed using the `filmpressure!` function
- `Fx :: Matrix{T}`: Total force acting on the fluid, x-component
- `Fy :: Matrix{T}`: Total force acting on the fluid, y-component
- `slipx :: Matrix{T}`: Friction force due to substrate slip, x-component
- `slipy :: Matrix{T}`: Friction force due to substrate slip, y-component
- `h∇px :: Matrix{T}`: Pressure gradient times the height, x-component
- `h∇py :: Matrix{T}`: Pressure gradient times the height, y-component
- `dgrad :: Array{T,N}`: Dummy allocation to store shifted arrays using `circshift!`

"""
Base.@kwdef struct State{T, N} <: LBM_state
    # Distribution functions
    fout :: Array{T, N}
    ftemp :: Array{T, N}
    feq :: Array{T, N} 
    # Macroscopic variables and moments 
    height :: Matrix{T} 
    velx :: Matrix{T} 
    vely :: Matrix{T}  
    vsq :: Matrix{T} 
    pressure :: Matrix{T} 
    # Forces and a dummy for the gradient
    Fx :: Matrix{T} 
    Fy :: Matrix{T} 
    slipx :: Matrix{T} 
    slipy :: Matrix{T} 
    h∇px :: Matrix{T} 
    h∇py :: Matrix{T} 
    dgrad :: Array{T,N} 
end

"""
    CuState

Data structure that stores all arrays for a given simulation.
Specific for GPU computing using `CUDA.jl`

# Arguments

- `fout :: Array{T,N}`: Output distribution function
- `ftemp :: Array{T,N}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: Array{T,N}`: Equilibrium distribution function
- `height :: Matrix{T}`: Field that stores the scalar height values
- `velx :: Matrix{T}`: Field that stores the x-component of the velocity vector
- `vely :: Matrix{T}`: Field that stores the y-component of the velocity vector
- `vsq :: Matrix{T}`: Field that stores the velocity squared, used in `equilibrium!`
- `pressure :: Matrix{T}`: Pressure distribution computed using the `filmpressure!` function
- `Fx :: Matrix{T}`: Total force acting on the fluid, x-component
- `Fy :: Matrix{T}`: Total force acting on the fluid, y-component
- `slipx :: Matrix{T}`: Friction force due to substrate slip, x-component
- `slipy :: Matrix{T}`: Friction force due to substrate slip, y-component
- `h∇px :: Matrix{T}`: Pressure gradient times the height, x-component
- `h∇py :: Matrix{T}`: Pressure gradient times the height, y-component
- `dgrad :: Array{T,N}`: Dummy allocation to store shifted arrays using `circshift!`

"""
struct CuState <: LBM_state
    # Distribution functions
    fout :: CuArray
    ftemp :: CuArray
    feq :: CuArray 
    # Macroscopic variables and moments 
    height :: CuArray
    velx :: CuArray
    vely :: CuArray
    vsq :: CuArray
    pressure :: CuArray
    # Forces and a dummy for the gradient
    Fx :: CuArray 
    Fy :: CuArray 
    slipx :: CuArray
    slipy :: CuArray 
    h∇px :: CuArray 
    h∇py :: CuArray 
    dgrad :: CuArray
end


"""
    StateActive{T, N}

Data structure that stores all arrays for a given simulation including a field rho that holds the distribution of active colloids in the fluid.

# Arguments

- `fout :: Array{T,N}`: Output distribution function
- `ftemp :: Array{T,N}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: Array{T,N}`: Equilibrium distribution function
- `height :: Matrix{T}`: Field that stores the scalar height values
- `velx :: Matrix{T}`: Field that stores the x-component of the velocity vector
- `vely :: Matrix{T}`: Field that stores the y-component of the velocity vector
- `vsq :: Matrix{T}`: Field that stores the velocity squared, used in `equilibrium!`
- `pressure :: Matrix{T}`: Pressure distribution computed using the `filmpressure!` function
- `Fx :: Matrix{T}`: Total force acting on the fluid, x-component
- `Fy :: Matrix{T}`: Total force acting on the fluid, y-component
- `slipx :: Matrix{T}`: Friction force due to substrate slip, x-component
- `slipy :: Matrix{T}`: Friction force due to substrate slip, y-component
- `h∇px :: Matrix{T}`: Pressure gradient times the height, x-component
- `h∇py :: Matrix{T}`: Pressure gradient times the height, y-component
- `dgrad :: Array{T,N}`: Dummy allocation to store shifted arrays using `circshift!`
- `rho :: Matrix{T}` : Distribution of active colloid density
- `lap_rho :: Matrix{T}` : Laplacian of rho
- `grad_rho_x :: Matrix{T}` : x compnent of gradient of rho
- `grad_rho_y :: Matrix{T}` : y compnent of gradient of rho
- `lap_h :: Matrix{T}` : Laplacian of height
- `grad_h_x :: Matrix{T}` : x compnent of gradient of height
- `grad_h_y :: Matrix{T}` : y compnent of gradient of height
- `rho_int :: Matrix{T}` : the derivative of rho
- `γ :: Matrix{T}` : surface tension field
- `θ :: Matrix{T}` : contact angle field
"""
Base.@kwdef struct StateActive{T, N} <: LBM_state
    # Distribution functions
    fout :: Array{T, N}
    ftemp :: Array{T, N}
    feq :: Array{T, N} 
    # Macroscopic variables and moments 
    height :: Matrix{T} 
    velx :: Matrix{T} 
    vely :: Matrix{T}  
    vsq :: Matrix{T} 
    pressure :: Matrix{T} 
    # Forces and a dummy for the gradient
    Fx :: Matrix{T} 
    Fy :: Matrix{T} 
    slipx :: Matrix{T} 
    slipy :: Matrix{T} 
    h∇px :: Matrix{T} 
    h∇py :: Matrix{T} 
    dgrad :: Array{T,N} 
    # active matter realted variables
    lap_rho :: Matrix{T}
    grad_rho_x :: Matrix{T}
    grad_rho_y :: Matrix{T}
    lap_h :: Matrix{T}
    grad_h_x :: Matrix{T}
    grad_h_y :: Matrix{T}
    rho :: Matrix{T}
    rho_int :: Matrix{T}
    γ :: Matrix{T}
    θ :: Matrix{T}
    k :: Matrix{T}
end

"""
    State_1D{T, N}

Data structure that stores all arrays for a given simulation.

# Arguments

- `fout :: Matrix{T}`: Output distribution function
- `ftemp :: Matrix{T}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: Matrix{T}`: Equilibrium distribution function
- `height :: Vector{T}`: Field that stores the scalar height values
- `vel :: Vector{T}`: Field that stores the velocity
- `pressure :: Vector{T}`: Pressure distribution computed using the `filmpressure!` function
- `F :: Vector{T}`: Total force acting on the fluid
- `slip :: Vector{T}`: Friction force due to substrate slip
- `h∇p :: Vector{T}`: Pressure gradient times the height
- `dgrad :: Matrix{T}`: Dummy allocation to store shifted arrays using `circshift!`
"""
Base.@kwdef struct State_1D{T} <: LBM_state_1D
    # Distribution functions
    fout :: Matrix{T}
    ftemp :: Matrix{T}
    feq :: Matrix{T} 
    # Macroscopic variables and moments 
    height :: Vector{T} 
    vel :: Vector{T} 
    pressure :: Vector{T} 
    # Forces and a dummy for the gradient
    F :: Vector{T} 
    slip :: Vector{T}
    h∇p :: Vector{T}
    dgrad :: Matrix{T} 
end


"""
    StateActive_1D{T, N}

Data structure that stores all arrays for a given simulation.

# Arguments

- `fout :: Matrix{T}`: Output distribution function
- `ftemp :: Matrix{T}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: Matrix{T}`: Equilibrium distribution function
- `height :: Vector{T}`: Field that stores the scalar height values
- `vel :: Vector{T}`: Field that stores the velocity
- `pressure :: Vector{T}`: Pressure distribution computed using the `filmpressure!` function
- `F :: Vector{T}`: Total force acting on the fluid
- `slip :: Vector{T}`: Friction force due to substrate slip
- `h∇p :: Vector{T}`: Pressure gradient times the height
- `dgrad :: Matrix{T}`: Dummy allocation to store shifted arrays using `circshift!`
- `rho :: Vector{T}` : Distribution of active colloid density
- `lap_rho :: Vector{T}` : Laplacian of rho
- `grad_rho :: Vector{T}` :  gradient of rho
- `lap_h :: Vector{T}` : Laplacian of height
- `grad_h :: Vector{T}` : x compnent of gradient of height
- `rho_int :: Vector{T}` : the derivative of rho
- `γ :: Vector{T}` : surface tension field
- `θ :: Vector{T}` : contact angle field
"""
Base.@kwdef struct StateActive_1D{T} <: Active_1D
    # Distribution functions
    fout :: Matrix{T}
    ftemp :: Matrix{T}
    feq :: Matrix{T} 
    gout :: Matrix{T}
    gtemp :: Matrix{T}
    geq :: Matrix{T} 
    gbound :: Matrix{T}
    hout :: Matrix{T}
    htemp :: Matrix{T}
    heq :: Matrix{T}
    bout :: Matrix{T}
    btemp :: Matrix{T}
    beq :: Matrix{T}
    # Macroscopic variables and moments 
    height :: Vector{T} 
    vel :: Vector{T} 
    pressure :: Vector{T} 
    γ :: Vector{T}
    ∇γ :: Vector{T}
    # Forces and a dummy for the gradient
    F :: Vector{T} 
    rho_F :: Vector{T} 
    slip :: Vector{T}
    h∇p :: Vector{T}
    dgrad :: Matrix{T} 
    # active matter realted variables
    rho :: Vector{T}
    grad_rho :: Vector{T}
    grad_h :: Vector{T}
    grad_p :: Vector{T}
    rho_vel :: Vector{T}
    rho_adv_vel :: Vector{T}
    rho_pressure :: Vector{T}
    z :: Vector{T}
    k :: Vector{T}
    stress :: Vector{T}
    ∇rho_u :: Vector{T}
    rho∇p :: Vector{T}
    slip_rho :: Vector{T}
    stress_rho :: Vector{T}
    diffusion_rho :: Vector{T}
    V_p :: Vector{T}
    V_gamma :: Vector{T}
    nabla_F :: Vector{T}
    rocks :: Vector{Int}
    rightrocks :: Vector{Int}
    leftrocks :: Vector{Int}
    D :: Vector{Float64}
    fourier_kernel :: Vector{Float64}
    # fields for the prduct concentration
    rho_A :: Vector{T}
    grad_rho_A :: Vector{T}
    rho_A_vel :: Vector{T}
    rho_A_adv_vel :: Vector{T}
    rho_A_pressure :: Vector{T}
    rho_A∇p :: Vector{T}
    rho_A_F :: Vector{T} 
    V_p_A :: Vector{T}
    V_gamma_A :: Vector{T}
    precursor :: Vector{T}
    rho_B :: Vector{T}
    grad_rho_B :: Vector{T}
    rho_B_vel :: Vector{T}
    rho_B_adv_vel :: Vector{T}
    rho_B_pressure :: Vector{T}
    rho_B∇p :: Vector{T}
    rho_B_F :: Vector{T} 
end



Base.@kwdef struct CuStateActive_1D{T} <: Active_1D
    # Distribution functions
    fout :: CuArray
    ftemp :: CuArray
    feq :: CuArray 
    gout :: CuArray
    gtemp :: CuArray
    geq :: CuArray 
    # Macroscopic variables and moments 
    height :: CuArray 
    vel :: CuArray 
    pressure :: CuArray 
    γ :: CuArray
    ∇γ :: CuArray
    # Forces and a dummy for the gradient
    F :: CuArray 
    rho_F :: CuArray
    slip :: CuArray
    h∇p :: CuArray
    # for gradients
    dgrad :: CuArray
    # active matter realted variables
    lap_rho :: CuArray
    grad_rho :: CuArray
    lap_h :: CuArray
    grad_h :: CuArray
    rho :: CuArray
    rho_vel :: CuArray
    rho_temp :: CuArray
    rho_int :: CuArray
    θ :: CuArray
    k :: CuArray
end


"""
    Sys(sysc, device, exotic)

Mostly allocations of arrays used to run a simulation, but all within one function :)

# Arguments

- `sysc :: SysConst`: Needed for the lattice dimensions, `Lx` and `Ly`
- `device :: String`: Use either `CPU` for computation of a CPU or `GPU` for computation on the GPU 
- `exotic :: Bool`: If true thermal fluctuations can be computed and saved to the `fthermalx` and `fthermaly` field
- `T <: Number`: Numerical type, it is strongly suggested to use `Float64`
"""
function Sys(sysc::SysConst, device::String, exotic::Bool, T)
    if device == "CPU"
        # Meso
        fout = zeros(sysc.Lx, sysc.Ly, 9)
        ftemp = zeros(sysc.Lx, sysc.Ly, 9)
        feq = zeros(sysc.Lx, sysc.Ly, 9)
        
        # Macro
        height = ones(sysc.Lx, sysc.Ly)
        velx = zeros(sysc.Lx, sysc.Ly)
        vely = zeros(sysc.Lx, sysc.Ly)
        vsq = zeros(sysc.Lx, sysc.Ly)
        pressure = zeros(sysc.Lx, sysc.Ly)
        dgrad = zeros(sysc.Lx, sysc.Ly, 8)

        # Forces
        Fx = zeros(sysc.Lx, sysc.Ly)
        Fy = zeros(sysc.Lx, sysc.Ly)
        slipx = zeros(sysc.Lx, sysc.Ly)
        slipy = zeros(sysc.Lx, sysc.Ly)
        h∇px = zeros(sysc.Lx, sysc.Ly)
        h∇py = zeros(sysc.Lx, sysc.Ly)
        if exotic
            # Probably many more forces here in the future
            fthermalx = zeros(sysc.Lx, sysc.Ly)
            fthermaly = zeros(sysc.Lx, sysc.Ly)
            return fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py, fthermalx, fthermaly
        end
        return fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py
    elseif device == "GPU"
        # Meso
        fout = CUDA.zeros(T, sysc.Lx, sysc.Ly, 9)
        ftemp = CUDA.zeros(T, sysc.Lx, sysc.Ly, 9)
        feq = CUDA.zeros(T, sysc.Lx, sysc.Ly, 9)
        
        # Macro
        height = CUDA.ones(T, sysc.Lx, sysc.Ly)
        velx = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        vely = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        vsq = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        pressure = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        dgrad = CUDA.zeros(T, sysc.Lx, sysc.Ly, 8)

        # Forces
        Fx = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        Fy = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        slipx = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        slipy = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        h∇px = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        h∇py = CUDA.zeros(T, sysc.Lx, sysc.Ly)
        if exotic
            fthermalx = CUDA.zeros(T, sysc.Lx, sysc.Ly)
            fthermaly = CUDA.zeros(T, sysc.Lx, sysc.Ly)
            return fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py, fthermalx, fthermaly
        end
        return fout, ftemp, feq, height, velx, vely, vsq, pressure, dgrad, Fx, Fy, slipx, slipy, h∇px, h∇py
    end
end

function Sys(sysc::SysConst, device::String; T=Float64)
    if device == "CPU"
        dyn = State{T, 3}(
            fout = zeros(sysc.Lx, sysc.Ly, 9),
            ftemp = zeros(sysc.Lx, sysc.Ly, 9),
            feq = zeros(sysc.Lx, sysc.Ly, 9),
            height = ones(sysc.Lx, sysc.Ly),
            velx = zeros(sysc.Lx, sysc.Ly),
            vely = zeros(sysc.Lx, sysc.Ly),
            vsq = zeros(sysc.Lx, sysc.Ly),
            pressure = zeros(sysc.Lx, sysc.Ly),
            dgrad = zeros(sysc.Lx, sysc.Ly, 8),
            Fx = zeros(sysc.Lx, sysc.Ly),
            Fy = zeros(sysc.Lx, sysc.Ly),
            slipx = zeros(sysc.Lx, sysc.Ly),
            slipy = zeros(sysc.Lx, sysc.Ly),
            h∇px = zeros(sysc.Lx, sysc.Ly),
            h∇py = zeros(sysc.Lx, sysc.Ly)
        )
        return dyn
    elseif device == "GPU"
        dyn = CuState(
            # Distribution functions
            CUDA.zeros(T, sysc.Lx, sysc.Ly, 9),
            CUDA.zeros(T, sysc.Lx, sysc.Ly, 9),
            CUDA.zeros(T, sysc.Lx, sysc.Ly, 9),
            # Macro
            CUDA.ones(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            # For gradients
            CUDA.zeros(T, sysc.Lx, sysc.Ly, 8),
            # Forces
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly),
            CUDA.zeros(T, sysc.Lx, sysc.Ly)
        )
        return dyn
    end
end
#initialising for active thin film model
"""
    Sys(sysc, device, exotic)

Mostly allocations of arrays used to run a simulation  for active thin film model

# Arguments

- `sysc :: SysConstActive`: Needed for the lattice dimensions, `Lx` and `Ly`
- `device :: String`: Use either `CPU` for computation of a CPU or `GPU` for computation on the GPU
- `T <: Number`: Numerical type, it is strongly suggested to use `Float64`
"""
function Sys(sysc::SysConstActive, device::String; T=Float64)    
    if device=="GPU"
        println("No CUDA available for active matter yet")
        # TODO: implement
    elseif device=="CPU"
        dyn = StateActive{T, 3}(
            fout = zeros(sysc.Lx, sysc.Ly, 9),
            ftemp = zeros(sysc.Lx, sysc.Ly, 9),
            feq = zeros(sysc.Lx, sysc.Ly, 9),
            height = ones(sysc.Lx, sysc.Ly),
            velx = zeros(sysc.Lx, sysc.Ly),
            vely = zeros(sysc.Lx, sysc.Ly),
            vsq = zeros(sysc.Lx, sysc.Ly),
            pressure = zeros(sysc.Lx, sysc.Ly),
            dgrad = zeros(sysc.Lx, sysc.Ly, 8),
            Fx = zeros(sysc.Lx, sysc.Ly),
            Fy = zeros(sysc.Lx, sysc.Ly),
            slipx = zeros(sysc.Lx, sysc.Ly),
            slipy = zeros(sysc.Lx, sysc.Ly),
            h∇px = zeros(sysc.Lx, sysc.Ly),
            h∇py = zeros(sysc.Lx, sysc.Ly),
            rho= ones(sysc.Lx, sysc.Ly),
            lap_rho = zeros(sysc.Lx, sysc.Ly),
            grad_rho_x = zeros(sysc.Lx, sysc.Ly),
            grad_rho_y = zeros(sysc.Lx, sysc.Ly),
            lap_h = zeros(sysc.Lx, sysc.Ly),
            grad_h_x = zeros(sysc.Lx, sysc.Ly),
            grad_h_y = zeros(sysc.Lx, sysc.Ly),
            rho_int = zeros(sysc.Lx, sysc.Ly),
            γ = ones(sysc.Lx, sysc.Ly) .* sysc.γ_0,
            θ = ones(sysc.Lx,sysc.Ly) .* sysc.θ_0,
            k = zeros(sysc.Lx, sysc.Ly)
        )
        return dyn
    end
end

"""
    Sys(sysc, exotic, T)

Mostly allocations of arrays used to run a simulation, but all within one function :)

# Arguments

- `sysc :: SysConst_1D`: Needed for the lattice dimensions, `L` 
- `exotic :: Bool`: If true thermal fluctuations can be computed and saved to the `fthermalx` and `fthermaly` field
- `T <: Number`: Numerical type, it is strongly suggested to use `Float64`
"""
function Sys(sysc::SysConst_1D, exotic::Bool, T)
    # Meso
    fout = zeros(sysc.L, 3)
    ftemp = zeros(sysc.L, 3)
    feq = zeros(sysc.L, 3)
    
    # Macro
    height = ones(sysc.L)
    vel = zeros(sysc.L)
    pressure = zeros(sysc.L)
    dgrad = zeros(sysc.L, 2)
    # Forces
    F = zeros(sysc.L)
    slip = zeros(sysc.L)
    h∇p = zeros(sysc.L)
    if exotic
        # Probably many more forces here in the future
        fthermal = zeros(sysc.L)
        
        return fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p, fthermal
    end

    return fout, ftemp, feq, height, vel, pressure, dgrad, F, slip, h∇p
    
end

function Sys(sysc::SysConst_1D; T=Float64)
    dyn = State_1D{T}(
        fout = zeros(sysc.L, 3),
        ftemp = zeros(sysc.L, 3),
        feq = zeros(sysc.L, 3),
        height = ones(sysc.L),
        vel = zeros(sysc.L),
        pressure = zeros(sysc.L),
        dgrad = zeros(sysc.L, 2),
        F = zeros(sysc.L),
        slip = zeros(sysc.L),
        h∇p = zeros(sysc.L) ,
    )
    return dyn   
end

#initialising active thin film model
"""
    Sys(sysc, device, exotic)

Mostly allocations of arrays used to run a simulation  for active thin film model

# Arguments

- `sysc :: SysConstActive_1D`: Needed for the lattice dimensions, `Lx` and `Ly`
- `device :: String`: Use either `CPU` for computation on a CPU or `GPU` for computation on the GPU
- `T <: Number`: Numerical type, it is strongly suggested to use `Float64`
"""
function Sys(sysc::SysConstActive_1D; T=Float64, device = "CPU")
    if device == "CPU"
        dyn = StateActive_1D{T}(
        fout = zeros(sysc.L, 3),
        ftemp = zeros(sysc.L, 3),
        feq = zeros(sysc.L, 3),
        gout = zeros(sysc.L, 3),
        gtemp = zeros(sysc.L, 3),
        geq = zeros(sysc.L, 3),
	    gbound = zeros(sysc.L, 3),
        hout = zeros(sysc.L, 3),
        heq = zeros(sysc.L, 3),
        htemp = zeros(sysc.L, 3),
        bout = zeros(sysc.L, 3),
        beq = zeros(sysc.L, 3),
        btemp = zeros(sysc.L, 3),
        height = ones(sysc.L),
        vel = zeros(sysc.L),
        rho_vel = zeros(sysc.L),
        pressure = zeros(sysc.L),
        dgrad = zeros(sysc.L, 2),
        F = zeros(sysc.L),
        slip = zeros(sysc.L),
        h∇p = zeros(sysc.L),
        rho= ones(sysc.L),
        grad_rho=zeros(sysc.L),
        grad_h=zeros(sysc.L),
        grad_p=zeros(sysc.L),
        γ = ones(sysc.L) .* sysc.γ_0,
        ∇γ = zeros(sysc.L),
        rho_F = zeros(sysc.L),
        rho_adv_vel=zeros(sysc.L),
        rho_pressure=zeros(sysc.L),
        z=zeros(sysc.L),
        k=zeros(sysc.L),
        stress=zeros(sysc.L),
        ∇rho_u=zeros(sysc.L),
        rho∇p=zeros(sysc.L),
        slip_rho=zeros(sysc.L),
        stress_rho=zeros(sysc.L),
        diffusion_rho=zeros(sysc.L),
	    V_p=zeros(sysc.L),
	    V_gamma=zeros(sysc.L),
	    nabla_F=zeros(sysc.L),
	    rocks=zeros(sysc.L),
	    rightrocks=zeros(sysc.L),
	    leftrocks=zeros(sysc.L),
        D= sysc.D .* ones(sysc.L),
        # fourier_kernel = [ sysc.production_rate /(sysc.D_Product * i^2 * ( 2 * pi / sysc.L )^(2) + sysc.sink_rate) for i in 0:Int(floor(sysc.L/2))]
        fourier_kernel = [ sysc.production_rate /(sysc.D_Product * i^2 * ( 2 * pi / sysc.L )^(2) + sysc.sink_rate_bulk) for i in 0:Int(floor(sysc.L/2))],
        rho_A = ones(sysc.L) ,
        grad_rho_A = zeros(sysc.L),
        rho_A_vel = zeros(sysc.L),
        rho_A_adv_vel = zeros(sysc.L),
        rho_A_pressure = zeros(sysc.L),
        rho_A∇p = zeros(sysc.L),
        rho_A_F = zeros(sysc.L),
        V_p_A=zeros(sysc.L),
	    V_gamma_A=zeros(sysc.L),
        precursor = zeros(sysc.L),
        rho_B = ones(sysc.L) ,
        grad_rho_B = zeros(sysc.L),
        rho_B_vel = zeros(sysc.L),
        rho_B_adv_vel = zeros(sysc.L),
        rho_B_pressure = zeros(sysc.L),
        rho_B∇p = zeros(sysc.L),
        rho_B_F = zeros(sysc.L),
        )
        return dyn
    elseif device=="GPU"
        dyn = CuStateActive_1D{T}(
             # Distribution functions
             CUDA.zeros(T, sysc.L, 3),
             CUDA.zeros(T, sysc.L, 3),
             CUDA.zeros(T, sysc.L, 3),
             CUDA.zeros(T, sysc.L, 3),
             CUDA.zeros(T, sysc.L, 3),
             CUDA.zeros(T, sysc.L, 3),
                # Macroscopic variables and moments 
             CUDA.ones(T, sysc.L),
             CUDA.zeros(T, sysc.L),
             CUDA.zeros(T, sysc.L),
             CUDA.zeros(T, sysc.L),
             CUDA.zeros(T, sysc.L),
                # Forces and a dummy for the gradient
            CUDA.zeros(T, sysc.L), 
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
                #Gradients
            CUDA.zeros(T,sysc.L,2),
                # active matter realted variables
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L),
            CUDA.zeros(T, sysc.L)
        )
        return dyn
    end
end
