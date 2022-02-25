# Parent type for all LBM states 
abstract type LBM_state end
# Derived type for 2D LBM states
abstract type LBM_state_2D <: LBM_state end

abstract type Expanded_2D <: LBM_state_2D end
# Derived type for 1D LBM states
abstract type LBM_state_1D <: LBM_state end

abstract type Expanded_1D <: LBM_state_1D end

# Parent type for system specific constants
abstract type Consts end

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
    SysConst{T}

Struct that contains the system size of a two dimensional system and the struct `Taumucs`.

# Arguments

- `Lx :: Int`: Length of the lattice sides in x-direction
- `Ly :: Int`: Length of the lattice sides in y-direction
- `s :: Taumucs` : Most of the run time constants

See also: [`SysConst_1D{T}`](@ref), [`Taumucs{T}`](@ref)

"""
Base.@kwdef struct SysConst{T} <: Consts
    # Lattice
    Lx :: Int = 256
    Ly :: Int = 256
    param :: Taumucs{T}
end

"""
    SysConst_1D{T}

Struct that contains the system size of a one dimensional system and the struct `Taumucs`.

# Arguments

- `L :: Int`: Number of lattice points
- `s :: Taumucs`: Most of the run time constants

See also: [`SysConst{T}`](@ref), [`Taumucs{T}`](@ref)

"""
Base.@kwdef struct SysConst_1D{T} <: Consts
    # Lattice
    L :: Int = 256
    param :: Taumucs{T}
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
Base.@kwdef struct State{T, N} <: LBM_state_2D
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
    State_thermal{T, N}

Data structure that contains a `State` and fields for thermal fluctuations.

# Arguments

- `basestate :: State{T,N}`: State data structer
- `kbtx :: Matrix{T}`: Force due to thermal fluctuations, x-component
- `kbty :: Matrix{T}`: Force due to thermal fluctuations, y-component

"""
Base.@kwdef struct State_thermal{T, N} <: Expanded_2D
    basestate :: State{T, N}
    kbtx :: Matrix{T} 
    kbty :: Matrix{T}  
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

struct CuState_thermal <: LBM_state_2D
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
    kbtx :: CuArray 
    kbty :: CuArray
    dgrad :: CuArray
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
    State_gamma_1D{T, N}

`State_1D` data structure with additional fields for surface tension `γ` and `∇γ`.

# Arguments

- `basestate :: State_1D{T}`: Base data structure for one dimensional simulations
- `γ :: Vector{T}`: Surface tension field
- `∇γ :: Vector{T}`: Surface tension gardient
"""
Base.@kwdef struct State_gamma_1D{T} <: Expanded_1D
    basestate :: State_1D{T}
    γ :: Vector{T}
    ∇γ :: Vector{T}
end

"""
    State_thermal_1D{T, N}

`State_1D` data structure with additional fields for thermal fluctuations (noise).

# Arguments

- `basestate :: State_1D{T}`: Base data structure for one dimensional simulations
- `kbt :: Vector{T}`: Surface tension field

"""
Base.@kwdef struct State_thermal_1D{T} <: Expanded_1D
    basestate :: State_1D{T}
    kbt :: Vector{T}
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
function Sys(sysc::SysConst, device::String; T=Float64, kind="simple")
    s = State{T, 3}(
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
                CUDA.zeros(T, sysc.Lx, sysc.Ly, 8)
            )
            return dyn
        end
    elseif kind == "thermal"
        if device == "CPU"
            dyn = State_thermal{T, 3}(
                basestate = s,
                kbtx = zeros(sysc.Lx, sysc.Ly),
                kbty = zeros(sysc.Lx, sysc.Ly)
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
                CUDA.zeros(T, sysc.Lx, sysc.Ly)
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
function Sys(sysc::SysConst_1D; T=Float64, kind="simple")
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
        dyn = State_thermal_1D{T}(
            basestate = s,
            kbt = zeros(sysc.L),
        )
    elseif kind == "gamma"
        dyn = State_gamma_1D{T}(
            basestate = s,
            γ = zeros(sysc.L),
            ∇γ = zeros(sysc.L)
        )
    end
    return dyn
    
end