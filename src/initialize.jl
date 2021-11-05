"""
    SysConst{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.
"""
@with_kw struct SysConst{T}
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
    g :: T = 0.0
end

"""
    SysConst_1D{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.
"""
@with_kw struct SysConst_1D{T}
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
    g :: T = 0.0
end

"""
    State{T, N}

Data structure that stores all arrays for a given simulation.
"""
@with_kw struct State{T, N}
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
"""
struct CuState
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
    State_1D{T, N}

Data structure that stores all arrays for a given simulation.
"""
@with_kw struct State_1D{T}
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
    Sys(sysc, device, exotic)

Mostly allocations of arrays used to run a simulation, but all within one function :)
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
        h∇p = zeros(sysc.L),
    )
    return dyn
    
end