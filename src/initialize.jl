# Parent type for all LBM states 
abstract type LBM_state end
# Derived type for 2D LBM states
abstract type LBM_state_2D <: LBM_state end

abstract type Expanded_2D <: LBM_state_2D end
abstract type MultiLayer_2D <: LBM_state_2D end
# Derived type for 1D LBM states
abstract type LBM_state_1D <: LBM_state end

abstract type Expanded_1D <: LBM_state_1D end
abstract type Boundary_1D <: Expanded_1D end
abstract type Active_1D <: LBM_state_1D end

# Parent type for system specific constants
abstract type Consts end

abstract type Consts_1D <: Consts end

"""
    Taumucs{T}

Struct that contains most run time constants, e.g. surface tension `γ`, viscosity `μ` and so on.

# Arguments

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

See also: [`SysConst{T}`](@ref), [`SysConst_1D{T}`](@ref)

"""
Base.@kwdef struct Taumucs{T} <: Consts
    Tmax::Int = 1000
    tdump::Int = Tmax ÷ 10
    # Collision related
    τ::T = 1.0
    cₛ::T = 1 / sqrt(3.0)
    μ::T = cₛ^2 * (τ - 0.5)
    # Force related
    δ::T = 1.0
    kbt::T = 0.0
    γ::T = 0.01
    n::Int = 9
    m::Int = 3
    hmin::T = 0.1
    hcrit::T = 0.05
    θ::T = 1 / 9
    g::T = 0.0
end

"""
    SysConst{T}

Struct that contains the system size of a two dimensional system and the struct `Taumucs`.

# Arguments

- `Lx :: Int`: Length of the lattice sides in x-direction
- `Ly :: Int`: Length of the lattice sides in y-direction
- `s :: Taumucs` : Most of the run time constants

See also: [`SysConst_1D{T}`](@ref), [`Taumucs{T}`](@ref)

"""
Base.@kwdef struct SysConst{T} <: Consts_1D
    # Lattice
    Lx::Int = 256
    Ly::Int = 256
    param::Taumucs{T}
end

"""
    SysConst_1D{T}

Struct that contains the system size of a one dimensional system and the struct `Taumucs`.

# Arguments

- `L :: Int`: Number of lattice points
- `s :: Taumucs`: Most of the run time constants

See also: [`SysConst{T}`](@ref), [`Taumucs{T}`](@ref)

"""
Base.@kwdef struct SysConst_1D{T} <: Consts_1D
    # Lattice
    L::Int = 256
    param::Taumucs{T}
end


"""
    SysConstWithBound_1D{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.

# Arguments


"""
Base.@kwdef struct SysConstWithBound_1D{T} <: Consts_1D
    # Lattice
    L::Int = 256
    param::Taumucs{T}
    obs::Array{T} = zeros(L)
    # interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld
    interior::Vector{T} = zeros(L)
    border = [zeros(L), zeros(L)]
end



"""
    SysConstMultiLayer{T}

Struct that contains all runtime constants for a multilayer lattice Boltzmann simulation. These constants define lattice parameters, physical properties, and numerical stabilizers used throughout the simulation. 

# Arguments

- `Lx :: Int`: Number of lattice points in the x-direction.
- `Ly :: Int`: Number of lattice points in the y-direction.
- `Tmax :: Int`: Total number of lattice Boltzmann time iterations.
- `tdump :: Int`: Interval for dumping data output during simulation, typically set as a fraction of `Tmax`.
- `layers :: Int`: Number of layers in the multilayer film model.

## Collision and Transport Properties
- `tau :: Array{T}`: Array of BGK relaxation rates, one for each layer.
- `cs :: T`: Lattice speed of sound, used as a reference velocity. Physical velocities should be less than `cs`.
- `mu :: Array{T}`: Array of kinematic viscosities for each layer, computed as `cs^2 * (tau - 0.5)`.
- `omeg :: Array{T}`: Collision kernel coefficients for interlayer interactions, with diagonal elements related to the relaxation rates.
- `umtau :: Array{T}`: Collision kernel weights, reciprocal of `tau` for self-interactions, zero for off-diagonal elements.

## Physical Forces and Interactions
- `delta :: Array{T}`: Slip lengths for each layer, defining how far the **no-slip** condition is interpolated into the substrate.
- `slip_coef :: Array{T}`: Coefficients for interlayer slip interactions
- `kbt :: T`: Thermal energy of the film; typically set to small values (e.g., 1e-7)
- `gamma :: Matrix{T}`: Surface tension coefficients, defined for layer interactions, We have gamma[i,j] the interaction between layer i and layer j where i=1 is the substrate, j=sys.layers+2 is the atmosphere, and we are only using the upper triangular matrix

## Disjoining Pressure and Stability
- `n :: Int`: Larger exponent used in the power-law representation of the disjoining pressure term.
- `m :: Int`: Smaller exponent used in the power-law representation of the disjoining pressure term.
- `hmin :: T`: Minimum film height at which the disjoining pressure functional vanishes. For simulations with three layers, a `hmin` of 0.75 or higher is advised to avoid instabilities in dewetting regions.
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term, used to regularize behavior near critical heights.

## Other Parameters
- `Tmax :: Int`: Total simulation time steps.
- `g :: T`: Gravitational acceleration, usually negligible in thin film simulations.
"""
Base.@kwdef struct SysConstMultiLayer{T} <: Consts
    # Lattice
    Lx :: Int = 256
    Ly :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    layers :: Int = 2
    # Collision related
    tau :: Array{T} = ones(layers)
    cs :: T = 1/sqrt(3.0)
    mu :: Array{T} = cs^2 .* (tau .- 0.5)
    omeg :: Array{T} = [i==j ? 1.0-1/tau[i] : 0.0 for i in 1:layers, j in 1:layers]
    umtau :: Array{T} = [i==j ? 1/tau[i] : 0.0 for i in 1:layers, j in 1:layers]
    # Force related
    delta :: Array{T} = ones(layers)
    # For test puruses
    slip_coef :: Array{T} = ones(layers, layers)
    kbt :: T = 0.0
    gamma :: Matrix{T} = 0.01 .* ones(layers+2, layers+2)
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
end


"""
    SysConstMultiLayer_1D{T}

Struct that contains all runtime constants for a 1D multilayer lattice Boltzmann simulation. These constants define lattice parameters, physical properties, and numerical stabilizers for the 1D simulation model.

# Arguments

- `L :: Int`: Number of lattice points.
- `Tmax :: Int`: Total number of lattice Boltzmann time iterations.
- `tdump :: Int`: Interval for dumping data output during simulation, typically set as a fraction of `Tmax`.
- `layers :: Int`: Number of layers in the multilayer film model.

## Collision and Transport Properties
- `tau :: Array{T}`: Array of BGK relaxation rates, one for each layer.
- `cs :: T`: Lattice speed of sound, used as a reference velocity. Physical velocities should be less than `cs`.
- `mu :: Array{T}`: Array of kinematic viscosities for each layer, computed as `cs^2 * (tau - 0.5)`.
- `omeg :: Array{T}`: Collision kernel coefficients for interlayer interactions, with diagonal elements related to the relaxation rates.
- `umtau :: Array{T}`: Collision kernel weights, reciprocal of `tau` for self-interactions, zero for off-diagonal elements.

## Physical Forces and Interactions
- `delta :: Array{T}`: Slip lengths for each layer, defining how far the **no-slip** condition is interpolated into the substrate.
- `slip_coef :: Array{T}`: Coefficients for interlayer slip interactions, primarily used for testing purposes.
- `kbt :: T`: Thermal energy of the film; typically set to small values (e.g., 1e-7).
- `gamma :: Matrix{T}`: Surface tension coefficients, defined for layer interactions.
- `repulsive :: Matrix{T}`: Repulsion factors to address the attraction that may arise at precursor height for positive spreading coefficients, based on `gamma`. It is one for negative spreading paramters and 1 for possitive spreading parameters

## Disjoining Pressure and Stability
- `n :: Int`: Larger exponent used in the power-law representation of the disjoining pressure term.
- `m :: Int`: Smaller exponent used in the power-law representation of the disjoining pressure term.
- `hmin :: T`: Minimum film height at which the disjoining pressure functional vanishes. 
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term, used to regularize behavior near critical heights.
- `extra_pressure_fac :: T`: Additional scaling factor for extra pressure contributions.

## Other Parameters
- `smooth_weight :: T`: Weight for a centered difference used to counteract errors that accumulate in the precursor layers over time.
- `usw :: T`: Smoothing weight factor defined as `1 - 2 * smooth_weight`.
"""

Base.@kwdef struct SysConstMultiLayer_1D{T} <: Consts
    # Lattice
    L :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    layers :: Int = 2
    # Collision related
    tau :: Array{T} = ones(layers)
    cs :: T = 1/sqrt(3.0)
    mu :: Array{T} = cs^2 .* (tau .- 0.5)
    omeg :: Array{T} = [i==j ? 1.0-1/tau[i] : 0.0 for i in 1:layers, j in 1:layers]
    umtau :: Array{T} = [i==j ? 1/tau[i] : 0.0 for i in 1:layers, j in 1:layers]
    # Force related
    delta :: Array{T} = ones(layers)
    # For test puruses
    slip_coef :: Array{T} = ones(layers, layers)
    kbt :: T = 0.0
    gamma :: Matrix{T} = 0.01 .* ones(layers+2, layers+2)
    repulsive :: Matrix{T} = repulsive_values(layers, gamma)
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    extra_pressure_fac :: T = 10.0
end

function repulsive_values(layers, gamma)
    rep=zeros(layers,layers)
    if layers==2
        rep[1,1]=(gamma[1,3]-gamma[1,2]-gamma[2,3])>0 ? -1.0 : 1.0
        rep[2,1]=(gamma[2,4]-gamma[2,3]-gamma[3,4])>0 ? -1.0 : 1.0
        rep[1,2]=(gamma[1,4]+gamma[2,3]-gamma[1,3]-gamma[2,4])>0 ? -1.0 : 1.0
    elseif layers==3
        rep[1,1]=(gamma[1,3]-gamma[1,2]-gamma[2,3])>0 ? 0.0 : 1.0
        rep[2,1]=(gamma[2,4]-gamma[2,3]-gamma[3,4])>0 ? 0.0 : 1.0
        rep[3,1]=(gamma[3,5]-gamma[3,4]-gamma[4,5])>0 ? 0.0 : 1.0
        # rep[1,2]=(gamma[1,4]+gamma[2,3]-gamma[1,3]-gamma[2,4])>0 ? 0.0 : 1.0
        # rep[2,2]=(gamma[2,5]+gamma[3,4]-gamma[2,4]-gamma[3,5])>0 ? 0.0 : 1.0
        # rep[1,3]=(gamma[1,5]+gamma[2,4]-gamma[1,4]-gamma[2,5])>0 ? 0.0 : 1.0
    else
        rep .= 1.0
    end
    return rep
end

"""
    SysConstMiscible_1D{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.

# Arguments

- `L :: Int`: Number of lattice points
- `Tmax :: Int`: Number of lattice Boltzmann time iterations
- `tdump :: Int`: Dumping interval for e.g. data output
- `liquids :: Int`: number of used liquids
- `tau :: T`: BGK relaxation rate 
- `umtau :: T`: 1/tau
- `cₛ :: T`: Lattice speed of sound, every physical velocity needs to be smaller than this! 
- `μ :: T`: Kinematic fluid viscosity
- `delta :: T`: Slip length, defines how far the **no-slip** condition is interpolated into the substrate
- `kbt :: T`: Thermal energy of the film, works with small values ≈ 10^(-7)
- `gamma :: Array{T}`: Surface tensions, liquid one: solid-vapour; gamma[3,1], solid-liquid; gamma[2,1], liquid-vapour; gamma[1,1], liquid two: solid-vapour; gamma[3,2], solid-liquid; gamma[2,2], solid-vapour; gamma[1,2]
- `n :: Int`: Greater exponent of the two used for the powerlaw in the disjoining pressure
- `m :: Int`: Smaller exponent of the two used for the powerlaw in the disjoining pressure
- `hmin :: T`: Height value at which the disjoining pressure functional vanishes
- `hcrit :: T`: Numerical stabilizer for the disjoining pressure term
- `θ :: T`: Contact angle in multiples of π
- `g :: T`: gravitational acceleration, usually neglected in thin film simulations
- `D :: T`: diffusion coefficient

"""
Base.@kwdef struct SysConstMiscible_1D{T} <: Consts
    # Lattice
    L :: Int = 256
    Tmax :: Int = 1000
    tdump :: Int = Tmax÷10
    liquids :: Int =2
    # Collision related
    tau :: T = 1.0
    umtau :: T = 1/tau
    omeg :: T = 1.0-1/tau
    cₛ :: T = 1/sqrt(3.0)
    μ :: T =  cₛ^2 *( tau - 0.5)
    # Force related
    delta :: T = 1.0
    kbt :: T = 0.0
    gamma :: Array{T} = 0.01 .* ones(3,liquids)
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    g :: T = 0.0
    D :: T = 1e-4
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
- `D :: T`: Diffusion constant of the catalyst
- `R :: T`: catalyst radius, not used
- `Γ :: T`: Surface tension change parameter product `rho_A`
- `rho_crit :: T`: parameter controling the precursor/disjoining pressure of rho see `rho_pressure!`
- `rho_n :: Int`: power of the precursor/disjoining pressure of rho see `rho_pressure!` deprecated
- `alpha :: T`: vertical distribution of catalyst, alpha<0 is concentrated at the liquid interface, alpha>0 is concentrated at the solid substrate, see [Richter et all](https://arxiv.org/abs/2402.14635)
- `A_2 :: T`: second order diffusion correction, not used very often
- `weight_center :: T`: moving averadge smoothing weight
- `production_rate :: T`: production rate
- `sigma_A_up :: T`: volatility rate of product
- `sigma_B_up :: T`: volatility rate of reactant
- `rho_BRes_sigma_B_down :: T`: source rate of reactant
- `D_Product :: T`: Diffusion coefficient of product, associated to `state.rho_A`
- `D_B :: T`:  Diffusion coefficient of reactant, associated to `state.rho_B`
- `rho_BRes :: T`: reservoir concentration of reactant
- `GammaB :: T`: surface tension change by reactant, `state.rho_B`
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
    rho_crit :: T = 0.05
    rho_n :: Int = 5
    alpha :: T = 0.0
    A_2 :: T = 0.0
    weight_center :: T = 1.0
    production_rate :: T =0.1
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

Data structure for both macroscopic variables and distribution functions.

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

See also: [`State_thermal{T, N}`](@ref), [`CuState`](@ref), [`State_1D{T, N}`](@ref)

"""
Base.@kwdef struct State{T,N} <: LBM_state_2D
    # Distribution functions
    fout::Array{T,N}
    ftemp::Array{T,N}
    feq::Array{T,N}
    # Macroscopic variables and moments 
    height::Matrix{T}
    velx::Matrix{T}
    vely::Matrix{T}
    vsq::Matrix{T}
    pressure::Matrix{T}
    # Forces and a dummy for the gradient
    Fx::Matrix{T}
    Fy::Matrix{T}
    slipx::Matrix{T}
    slipy::Matrix{T}
    h∇px::Matrix{T}
    h∇py::Matrix{T}
    dgrad::Array{T,N}
end

"""
    State_thermal{T, N}

Data structure that contains a `State` and fields for thermal fluctuations.

# Arguments

- `basestate :: State{T,N}`: State data structer
- `kbtx :: Matrix{T}`: Force due to thermal fluctuations, x-component
- `kbty :: Matrix{T}`: Force due to thermal fluctuations, y-component

"""
Base.@kwdef struct State_thermal{T,N} <: Expanded_2D
    basestate::State{T,N}
    kbtx::Matrix{T}
    kbty::Matrix{T}
end

"""
    CuState

`State` data structure but for CUDA like memory.

Similar to `CuState_thermal` without the thermal fields.

# Arguments

- `fout :: CuArray{T,N}`: Output distribution function
- `ftemp :: CuArray{T,N}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: CuArray{T,N}`: Equilibrium distribution function
- `height :: CuArray{T}`: Field that stores the scalar height values
- `velx :: CuMatrix{T}`: Field that stores the x-component of the velocity vector
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
struct CuState <: LBM_state_2D
    # Distribution functions
    fout::CuArray
    ftemp::CuArray
    feq::CuArray
    # Macroscopic variables and moments 
    height::CuArray
    velx::CuArray
    vely::CuArray
    vsq::CuArray
    pressure::CuArray
    # Forces and a dummy for the gradient
    Fx::CuArray
    Fy::CuArray
    slipx::CuArray
    slipy::CuArray
    h∇px::CuArray
    h∇py::CuArray
    dgrad::CuArray
end

struct CuState_thermal <: LBM_state_2D
    # Distribution functions
    fout::CuArray
    ftemp::CuArray
    feq::CuArray
    # Macroscopic variables and moments 
    height::CuArray
    velx::CuArray
    vely::CuArray
    vsq::CuArray
    pressure::CuArray
    # Forces and a dummy for the gradient
    Fx::CuArray
    Fy::CuArray
    slipx::CuArray
    slipy::CuArray
    h∇px::CuArray
    h∇py::CuArray
    kbtx::CuArray
    kbty::CuArray
    dgrad::CuArray
end

"""
    State_1D{T, N}

Data structure that for a one dimensional simulation.

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
    fout::Matrix{T}
    ftemp::Matrix{T}
    feq::Matrix{T}
    # Macroscopic variables and moments 
    height::Vector{T}
    vel::Vector{T}
    pressure::Vector{T}
    # Forces and a dummy for the gradient
    F::Vector{T}
    slip::Vector{T}
    h∇p::Vector{T}
    dgrad::Matrix{T}
end


"""
    State_gamma_1D{T, N}

`State_1D` data structure with additional fields for surface tension `γ` and `∇γ`.

# Arguments

- `basestate :: State_1D{T}`: Base data structure for one dimensional simulations
- `γ :: Vector{T}`: Surface tension field
- `∇γ :: Vector{T}`: Surface tension gardient
"""
Base.@kwdef struct State_gamma_1D{T} <: Expanded_1D
    basestate::State_1D{T}
    γ::Vector{T}
    ∇γ::Vector{T}
end


"""
    StateWithBound_1D{T, N}

Data structure that stores all arrays for a given simulation.

"""
Base.@kwdef struct StateWithBound_1D{T} <: Boundary_1D
    basestate::State_1D{T}
    # surface tension gradient
    γ::Vector{T}
    ∇γ::Vector{T}
    # Distribution functions 
    fbound::Matrix{T}
end


"""
    State_thermal_1D{T, N}

`State_1D` data structure with additional fields for thermal fluctuations (noise).

# Arguments

- `basestate :: State_1D{T}`: Base data structure for one dimensional simulations
- `kbt :: Vector{T}`: Surface tension field

"""
Base.@kwdef struct State_thermal_1D{T} <: Expanded_1D
    basestate::State_1D{T}
    kbt::Vector{T}
end

"""
    CuStateMultiLayer{T}

Data structure for storing all arrays required in a multilayer lattice Boltzmann simulation, designed for GPU computation using `CuArray`. This structure includes distribution functions, macroscopic variables, forces, and auxiliary fields utilized during the simulation.

# Arguments

## Distribution Functions
- `fout :: CuArray{T}`: Output distribution function.
- `ftemp :: CuArray{T}`: Temporary distribution function, used if `sys.τ ≠ 1`.
- `feq :: CuArray{T}`: Equilibrium distribution function.

## Macroscopic Variables
- `height :: CuArray{T}`: Scalar height values for the fluid layers.
- `velx :: CuArray{T}`: Velocity field in the x-direction.
- `vely :: CuArray{T}`: Velocity field in the y-direction.
- `vsq :: CuArray{T}`: Velocity squared, used for computing equilibrium distributions.
- `pressure :: CuArray{T}`: Pressure distribution, computed using the `filmpressure!` function.

## Forces and Friction
- `Fx :: CuArray{T}`: Total force acting on the fluid in the x-direction.
- `Fy :: CuArray{T}`: Total force acting on the fluid in the y-direction.
- `slipx :: CuArray{T}`: Friction force due to substrate slip in the x-direction.
- `slipy :: CuArray{T}`: Friction force due to substrate slip in the y-direction.

## Gradients and Derived Quantities
- `h∇px :: CuArray{T}`: Pressure gradient multiplied by the height, in the x-direction.
- `h∇py :: CuArray{T}`: Pressure gradient multiplied by the height, in the y-direction.
- `dgrad :: CuArray{T}`: Dummy allocation for storing shifted arrays using `circshift!`.
- `hi :: CuArray{T}`: Sum of heights across layers (e.g., `h_1 + h_2`), used in multilayer calculations.
- `grad_hx :: CuArray{T}`: Gradient of the height field in the x-direction.
- `grad_hy :: CuArray{T}`: Gradient of the height field in the y-direction.
- `grad_h_sq :: CuArray{T}`: Squared gradient of the height field, used in pressure calculations.
"""
Base.@kwdef struct CuStateMultiLayer{T} <: MultiLayer_2D
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
    hi :: CuArray
    grad_hx :: CuArray
    grad_hy :: CuArray
    grad_h_sq :: CuArray
end

"""
    StateMultiLayer_1D{T}

Data structure for storing all arrays required for a 1D multilayer lattice Boltzmann simulation. This structure includes distribution functions, macroscopic variables, forces, and auxiliary fields utilized during the simulation.

# Arguments

## Distribution Functions
- `fout :: Matrix{T}`: Output distribution function for the simulation.
- `ftemp :: Matrix{T}`: Temporary distribution function, used if `sys.τ ≠ 1`.
- `feq :: Matrix{T}`: Equilibrium distribution function.

## Macroscopic Variables
- `height :: Vector{T}`: Scalar height values for the fluid layers.
- `vel :: Vector{T}`: Velocity field in 1D.
- `pressure :: Vector{T}`: Pressure distribution, computed using the `filmpressure!` function.

## Forces and Friction
- `F :: Vector{T}`: Total force acting on the fluid in 1D.
- `slip :: Vector{T}`: Friction force due to substrate slip.
- `h∇p :: Vector{T}`: Pressure gradient multiplied by the height.

## Gradients and Auxiliary Fields
- `dgrad :: Matrix{T}`: Dummy allocation for storing shifted arrays using `circshift!`.
- `hi :: Vector{T}`: Sum of heights across layers (e.g., `h_1 + h_2`), used in multilayer calculations.
- `grad_h :: Vector{T}`: Gradient of the height field, used in pressure and force calculations.
"""

Base.@kwdef struct StateMultiLayer_1D{T} <: LBM_state_1D
    # Distribution functions
    fout :: Array{T}
    ftemp :: Array{T}
    feq :: Array{T} 
    # Macroscopic variables and moments 
    height :: Array{T} 
    vel :: Array{T} 
    pressure :: Array{T} 
    # Forces and a dummy for the gradient
    F :: Array{T} 
    slip :: Array{T}
    h∇p :: Array{T}
    dgrad :: Array{T} 
    hi :: Array{T}
    grad_h :: Array{T}
end


"""
    StateActive_1D{T}

Data structure that stores all arrays for a given simulation of active 1D systems, including distribution functions, macroscopic variables, forces, and various fields associated with the system's behavior.

# Arguments

- `fout :: Matrix{T}`: Output distribution function for the film height `state.height`
- `ftemp :: Matrix{T}`: Temporary distribution function for the film height `state.height`
- `feq :: Matrix{T}`: Equilibrium distribution function for the film height `state.height`
- `gout :: Matrix{T}`: Output distribution function for the catalyst concentration `state.rho`
- `gtemp :: Matrix{T}`: Temporary distribution function for the catalyst concentration `state.rho`
- `geq :: Matrix{T}`: Equilibrium distribution function for the catalyst concentration `state.rho`
- `hout :: Matrix{T}`: Output distribution function for the product concentration `state.rho_A`
- `htemp :: Matrix{T}`: Temporary distribution function for the product concentration `state.rho_A`
- `heq :: Matrix{T}`: Equilibrium distribution function for the product concentration `state.rho_A`
- `bout :: Matrix{T}`: Output distribution function for the reactant concentration `state.rho_B`
- `btemp :: Matrix{T}`: Temporary distribution function for the reactant concentration `state.rho_B`
- `beq :: Matrix{T}`: Equilibrium distribution function for the reactant concentration `state.rho_B`


- `height :: Vector{T}`: Field that stores the scalar height values of the system
- `vel :: Vector{T}`: Field that stores the velocity of the system
- `pressure :: Vector{T}`: Pressure distribution, computed using the `filmpressure!` function
- `γ :: Vector{T}`: Surface tension field
- `∇γ :: Vector{T}`: Gradient of the surface tension field

- `F :: Vector{T}`: Total force acting on the fluid
- `rho_F :: Vector{T}`: Distribution of force density (related to colloid density and force)
- `slip :: Vector{T}`: Friction force due to substrate slip
- `h∇p :: Vector{T}`: Pressure gradient times the height

- `dgrad :: Matrix{T}`: Dummy allocation to store shifted arrays using `circshift!`
- `rho :: Vector{T}`: Distribution of catalyst density
- `grad_rho :: Vector{T}`: Gradient of the active colloid density field
- `grad_h :: Vector{T}`: Gradient of the height field (in the x-direction)
- `grad_p :: Vector{T}`: Gradient of the pressure field
- `rho_vel :: Vector{T}`: Velocity of the catalyst
- `rho_adv_vel :: Vector{T}`: Advected velocity field for catalyst
- `rho_pressure :: Vector{T}`: Pressure distribution for the catalyst phase

- `z :: Vector{T}`: Auxiliary variable (possibly related to concentration or other system variables)
- `k :: Vector{T}`: Constant related to the active matter system
- `stress :: Vector{T}`: Stress field in the system, derived from the surface tension gradient nabla gamm, see [Richter et all](https://arxiv.org/abs/2402.14635)
- `∇rho_u :: Vector{T}`: Gradient of catalyst density multiplied by the velocity field
- `rho∇p :: Vector{T}`: Product of catalyst density and pressure gradient
- `slip_rho :: Vector{T}`: Friction force due to slip and catalyst density
- `stress_rho :: Vector{T}`: Stress field related to catalyst density
- `diffusion_rho :: Vector{T}`: Diffusion term for catalyst density
- `V_p :: Vector{T}`: Potential field related to pressure
- `V_gamma :: Vector{T}`: Potential field related to surface tension
- `nabla_F :: Vector{T}`: Gradient of the total force field

- `D :: Vector{Float64}`: Diffusion coefficients catalyst
- `fourier_kernel :: Vector{Float64}`: Fourier transform kernel used in system analysis or simulations

# Product Concentration Fields

- `rho_A :: Vector{T}`: Distribution of active catalyst density for product `state.rho_A`
- `grad_rho_A :: Vector{T}`: Gradient of the density field for product `state.rho_A`
- `rho_A_vel :: Vector{T}`: Velocity field for product `state.rho_A`
- `rho_A_adv_vel :: Vector{T}`: Advected velocity field for product `state.rho_A`
- `rho_A_pressure :: Vector{T}`: Pressure distribution for product `state.rho_A`
- `rho_A∇p :: Vector{T}`: Product of catalyst density for product `state.rho_A` and pressure gradient
- `rho_A_F :: Vector{T}`: Total force on product `state.rho_A`
- `V_p_A :: Vector{T}`: Potential related to pressure for product `state.rho_A`
- `V_gamma_A :: Vector{T}`: Potential related to surface tension for product `state.rho_A`
- `precursor :: Vector{T}`: Precursor field for product `state.rho_A`, telling you where the film is dewetted and where not. 

- `rho_B :: Vector{T}`: Distribution of active catalyst density for reactant `state.rho_B`
- `grad_rho_B :: Vector{T}`: Gradient of the density field for reactant `state.rho_B`
- `rho_B_vel :: Vector{T}`: Velocity field for reactant `state.rho_B`
- `rho_B_adv_vel :: Vector{T}`: Advected velocity field for reactant `state.rho_B`
- `rho_B_pressure :: Vector{T}`: Pressure distribution for reactant `state.rho_B`
- `rho_B∇p :: Vector{T}`: Product of colloid density for reactant `state.rho_B` and pressure gradient
- `rho_B_F :: Vector{T}`: Total force on reactant `state.rho_B`
"""
Base.@kwdef struct StateActive_1D{T} <: Active_1D
    # Distribution functions
    fout :: Matrix{T}
    ftemp :: Matrix{T}
    feq :: Matrix{T}
    gout :: Matrix{T}
    gtemp :: Matrix{T}
    geq :: Matrix{T}
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
    
    # Active matter related variables
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
    
    # Diffusion and Fourier kernel fields
    D :: Vector{Float64}
    fourier_kernel :: Vector{Float64}
    
    # Fields for product concentration
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

"""
    StateMiscible_1D{T, N}

Data structure that stores all arrays for a given simulation.

# Arguments

- `fout :: Matrix{T}`: Output distribution function, 3x2 array
- `ftemp :: Matrix{T}`: Temporary distribution function, 3x2 array
- `feq :: Matrix{T}`: Equilibrium distribution function, 3x2 array
- `height :: Vector{T}`: Field that stores the scalar height values
- `vel :: Vector{T}`: Field that stores the velocity
- `pressure :: Vector{T}`: Pressure distribution computed using the `filmpressure!` function
- `F :: Vector{T}`: Total force acting on the fluid
- `slip :: Vector{T}`: Friction force due to substrate slip
- `h∇p :: Vector{T}`: Pressure gradient times the height
- `dgrad :: Matrix{T}`: Dummy allocation to store shifted arrays using `circshift!`
- `gamma :: Matrix{T}`: local surface tension field
- `stress :: Matrix{T}`: Marangoni shear stress
- `∇gamma :: Matrix{T}`: gradient gamma
- `∇grad_h :: Matrix{T}`: gradient gamma
"""
Base.@kwdef struct StateMiscible_1D{T} <: LBM_state_1D
    # Distribution functions
    fout :: Array{T}
    ftemp :: Array{T}
    feq :: Array{T} 
    # Macroscopic variables and moments 
    height :: Array{T} 
    vel :: Array{T} 
    pressure :: Array{T} 
    # Forces and a dummy for the gradient
    F :: Array{T} 
    slip :: Array{T}
    h∇p :: Array{T}
    dgrad :: Array{T} 
    gamma :: Array{T}
    stress :: Array{T}
    ∇gamma :: Vector{T}
    grad_h :: Array{T}
end

"""
    State_curved_1D{T, N}

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
- `substrate :: Vector{T}`: substrate profile
- `grad_substrate :: Vector{T}`: gradient of the substrate profile, needed for the pressure
"""
Base.@kwdef struct State_curved_1D{T} <: LBM_state_1D
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
    #substrate field
    substrate :: Vector{T}
    grad_substrate :: Vector{T}
end





"""
    Sys(sysc, device, exotic, T)

Mostly allocations of arrays used to run a simulation, but all within one function :)

Returns not a data structure such as state, but every array.
Therefore it is somewhat outdated by now (2022).

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
            return fout,
            ftemp,
            feq,
            height,
            velx,
            vely,
            vsq,
            pressure,
            dgrad,
            Fx,
            Fy,
            slipx,
            slipy,
            h∇px,
            h∇py,
            fthermalx,
            fthermaly
        end
        return fout,
        ftemp,
        feq,
        height,
        velx,
        vely,
        vsq,
        pressure,
        dgrad,
        Fx,
        Fy,
        slipx,
        slipy,
        h∇px,
        h∇py
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
            return fout,
            ftemp,
            feq,
            height,
            velx,
            vely,
            vsq,
            pressure,
            dgrad,
            Fx,
            Fy,
            slipx,
            slipy,
            h∇px,
            h∇py,
            fthermalx,
            fthermaly
        end
        return fout,
        ftemp,
        feq,
        height,
        velx,
        vely,
        vsq,
        pressure,
        dgrad,
        Fx,
        Fy,
        slipx,
        slipy,
        h∇px,
        h∇py
    end
end

"""
    Sys(sysc, device; T, kind)

Allocations of arrays used to populate the `State` data structure.

Returns a `State` data structure based on `sysc`, either one dimensional or two dimensional.

# Arguments

- `sysc :: SysConst`: Needed for the lattice dimensions, `Lx` and `Ly`
- `device :: String`: Use either `CPU` for computation of a CPU or `GPU` for computation on the GPU (GPU can only be used with a two dimensional system)
- `T <: Number`: Numerical type, it is strongly suggested to use `Float64`
- `kind :: String`: Indicator for different `State`'s, default value is "simple" which creates a `State` data structure, valid options are ["simple", "thermal"] 
"""
function Sys(sysc::SysConst, device::String; T = Float64, kind = "simple")
    s = State{T,3}(
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
    )
    if kind == "simple"
        if device == "CPU"
            dyn = s
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
                # Forces
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                # For gradients
                CUDA.zeros(T, sysc.Lx, sysc.Ly, 8),
            )
            return dyn
        end
    elseif kind == "thermal"
        if device == "CPU"
            dyn = State_thermal{T,3}(
                basestate = s,
                kbtx = zeros(sysc.Lx, sysc.Ly),
                kbty = zeros(sysc.Lx, sysc.Ly),
            )
            return dyn
        elseif device == "GPU"
            dyn = CuState_thermal(
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
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
                CUDA.zeros(T, sysc.Lx, sysc.Ly),
            )
            return dyn
        end
    end
end

"""
    Sys(sysc, T, kind)

Allocations of arrays used to run a one dimensional simulation.

Returns a `State` data structure based on `kind` the struct can be "simple", "thermal" or "gamma".

# Arguments

- `sysc :: SysConst_1D`: Stores the lattice size `L` 
- `T <: Number`: Optional numerical type, default is set to `Float64`
- `kind :: String`: Optional, default is set to `simple`
"""
function Sys(sysc::Consts_1D; T = Float64, kind = "simple")
    s = State_1D{T}(
        fout = zeros(sysc.L, 3),
        ftemp = zeros(sysc.L, 3),
        feq = zeros(sysc.L, 3),
        height = ones(sysc.L),
        vel = zeros(sysc.L),
        pressure = zeros(sysc.L),
        dgrad = zeros(sysc.L, 2),
        F = zeros(sysc.L),
        slip = zeros(sysc.L),
        h∇p = zeros(sysc.L),
    )
    if kind == "simple"
        dyn = s
    elseif kind == "thermal"
        dyn = State_thermal_1D{T}(basestate = s, kbt = zeros(sysc.L))
    elseif kind == "gamma"
        dyn = State_gamma_1D{T}(basestate = s, γ = zeros(sysc.L), ∇γ = zeros(sysc.L))
    elseif kind == "gamma_bound"
        dyn = StateWithBound_1D{T}(
            basestate = s,
            γ = zeros(sysc.L),
            ∇γ = zeros(sysc.L),
            fbound = zeros(sysc.L, 3),
        )
    elseif kind == "curved"
            if device=="CPU"
            	    dyn = State_curved_1D{T}(
                	    fout = zeros(sysc.L, 3),
                	    ftemp = zeros(sysc.L, 3),
                	    feq = zeros(sysc.L, 3),
                	    height = ones(sysc.L),
                	    vel = zeros(sysc.L),
                	    pressure = zeros(sysc.L),
                	    dgrad = zeros(sysc.L, 2),
               		    F = zeros(sysc.L),
                	    slip = zeros(sysc.L),
                	    h∇p = zeros(sysc.L),
                	    substrate = zeros(sysc.L),
                	    grad_substrate = zeros(sysc.L)
            	    )
            	    return dyn
            else
            	    println("Device not available")
            end
    end
    return dyn

end


function Sys(sysc::SysConstActive_1D; T=Float64, device = "CPU")
    if device == "CPU"
        dyn = StateActive_1D{T}(
        fout = zeros(sysc.L, 3),
        ftemp = zeros(sysc.L, 3),
        feq = zeros(sysc.L, 3),
        gout = zeros(sysc.L, 3),
        gtemp = zeros(sysc.L, 3),
        geq = zeros(sysc.L, 3),
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
        D= sysc.D .* ones(sysc.L),
        fourier_kernel = ones(sysc.L),
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
	    println("This device is not available")
        return dyn
    end
end


function Sys(sysc::SysConstMultiLayer; T=Float64, device="CPU", kind="")
    if kind ==""
        if device=="CPU"
            dyn = StateMultiLayer{T}(
                fout = zeros(sysc.Lx, sysc.Ly, 9, sysc.layers),
                ftemp = zeros(sysc.Lx, sysc.Ly, 9, sysc.layers),
                feq = zeros(sysc.Lx, sysc.Ly, 9, sysc.layers),
                height = ones(sysc.Lx, sysc.Ly, sysc.layers),
                velx = zeros(sysc.Lx, sysc.Ly, sysc.layers),
                vely = zeros(sysc.Lx, sysc.Ly, sysc.layers),
                vsq = zeros(sysc.Lx, sysc.Ly, sysc.layers),
                pressure = zeros(sysc.Lx, sysc.Ly, sysc.layers),
                dgrad = zeros(sysc.Lx, sysc.Ly, sysc.layers, 8),
                Fx = zeros(sysc.Lx, sysc.Ly, sysc.layers),
                Fy = zeros(sysc.Lx, sysc.Ly, sysc.layers),
                slipx = zeros(sysc.Lx,sysc.Ly, sysc.layers),
                slipy = zeros(sysc.Lx,sysc.Ly, sysc.layers),
                h∇px = zeros(sysc.Lx,sysc.Ly, sysc.layers),
                h∇py = zeros(sysc.Lx,sysc.Ly, sysc.layers),
                hi = zeros(sysc.Lx,sysc.Ly, sysc.layers -1),
                grad_hx = zeros(sysc.Lx,sysc.Ly, sysc.layers -1),
                grad_hy = zeros(sysc.Lx,sysc.Ly, sysc.layers -1),
                grad_h_sq = zeros(sysc.Lx,sysc.Ly, sysc.layers -1)
            )
            return dyn
        elseif device=="GPU"
            dyn = CuStateMultiLayer{T}(
                fout = CUDA.zeros(T,sysc.Lx, sysc.Ly, 9, sysc.layers),
                ftemp = CUDA.zeros(T,sysc.Lx, sysc.Ly, 9, sysc.layers),
                feq = CUDA.zeros(T,sysc.Lx, sysc.Ly, 9, sysc.layers),
                height = CUDA.ones(T,sysc.Lx, sysc.Ly, sysc.layers),
                velx = CUDA.zeros(T,sysc.Lx, sysc.Ly, sysc.layers),
                vely = CUDA.zeros(T,sysc.Lx, sysc.Ly, sysc.layers),
                vsq = CUDA.zeros(T,sysc.Lx, sysc.Ly, sysc.layers),
                pressure = CUDA.zeros(T,sysc.Lx, sysc.Ly, sysc.layers),
                dgrad = CUDA.zeros(T,sysc.Lx, sysc.Ly,sysc.layers,8),
                Fx = CUDA.zeros(T,sysc.Lx, sysc.Ly, sysc.layers),
                Fy = CUDA.zeros(T,sysc.Lx, sysc.Ly, sysc.layers),
                slipx = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers),
                slipy = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers),
                h∇px = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers),
                h∇py = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers),
                hi = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers -1),
                grad_hx = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers -1),
                grad_hy = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers -1),
                grad_h_sq = CUDA.zeros(T,sysc.Lx,sysc.Ly, sysc.layers -1)
            )
             return dyn
        else
            throw(DomainError(device, "Unknown device \"$device\". returning nothing"))
        end
    end
end




function Sys(sysc::SysConstMiscible_1D; T=Float64, device="CPU", kind="")
    if kind ==""
        if device=="CPU"
            dyn = StateMiscible_1D{T}(
                fout = zeros(sysc.L, 3, sysc.liquids),
                ftemp = zeros(sysc.L, 3, sysc.liquids),
                feq = zeros(sysc.L, 3, sysc.liquids),
                height = ones(sysc.L, sysc.liquids),
                vel = zeros(sysc.L, sysc.liquids),
                pressure = zeros(sysc.L, sysc.liquids),
                dgrad = zeros(sysc.L, 2, sysc.liquids),
                F = zeros(sysc.L, sysc.liquids),
                slip = zeros(sysc.L, sysc.liquids),
                h∇p = zeros(sysc.L, sysc.liquids),
                gamma = zeros(sysc.L, 3),
                stress = zeros(sysc.L, sysc.liquids),
                ∇gamma = zeros(sysc.L),
                grad_h = zeros(sysc.L, sysc.liquids)
            )
            return dyn
        end
    end
end



function Sys(sysc::SysConstMultiLayer_1D; T=Float64, device="CPU", kind="")
    if kind ==""
        if device=="CPU"
            dyn = StateMultiLayer_1D{T}(
                fout = zeros(sysc.L, 3, sysc.layers),
                ftemp = zeros(sysc.L, 3, sysc.layers),
                feq = zeros(sysc.L, 3, sysc.layers),
                height = ones(sysc.L, sysc.layers),
                vel = zeros(sysc.L, sysc.layers),
                pressure = zeros(sysc.L, sysc.layers),
                dgrad = zeros(sysc.L, 2, sysc.layers),
                F = zeros(sysc.L, sysc.layers),
                slip = zeros(sysc.L, sysc.layers),
                h∇p = zeros(sysc.L, sysc.layers),
                hi = zeros(sysc.L, sysc.layers),
                grad_h = zeros(sysc.L, sysc.layers -1)
            )
            return dyn
        end
    end
end


