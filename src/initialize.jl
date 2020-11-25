"""
    SysConst{T}

Struct that contains all run time constants, e.g. lattice size, surface tension `γ` and so on.
"""
@with_kw struct SysConst{T}
    # Lattice
    Lx :: Int = 256
    Ly :: Int = 256
    Tmax :: Int = 1000
    # Collision related
    τ :: T = 1.0
    cₛ :: T = 1/sqrt(3.0)
    μ :: T =  cₛ^2 / 2 *( τ - 0.5)
    # Force related
    kbt :: T = 0.0
    γ :: T = 0.01
    n :: Int = 9
    m :: Int = 3
    hmin :: T = 0.1
    hcrit :: T = 0.05
    g :: T = 0.0
end

"""
    Sys(sysc, device, exotic)

Mostly allocations of arrays used to run a simulation, but all within one function :)
"""
function Sys(sysc :: SysConst, device::String, exotic::Bool)
    if device == "CPU"
        # Meso
        fout = zeros(sysc.Lx, sysc.Ly, 9)
        ftemp = zeros(sysc.Lx, sysc.Ly, 9)
        feq = zeros(sysc.Lx, sysc.Ly, 9)
        
        # Macro
        height = ones(sysc.Lx, sysc.Ly)
        velx = zeros(sysc.Lx, sysc.Ly)
        vely = zeros(sysc.Lx, sysc.Ly)
        pressure = zeros(sysc.Lx, sysc.Ly)

        # Forces
        Fx = zeros(sysc.Lx, sysc.Ly)
        Fy = zeros(sysc.Lx, sysc.Ly)
        slipx = zeros(sysc.Lx, sysc.Ly)
        slipy = zeros(sysc.Lx, sysc.Ly)
        h∇px = zeros(sysc.Lx, sysc.Ly)
        h∇py = zeros(sysc.Lx, sysc.Ly)
        if exotic
            fthermalx = zeros(sysc.Lx, sysc.Ly)
            fthermaly = zeros(sysc.Lx, sysc.Ly)
            return fout, ftemp, feq, height, velx, vely, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py, fthermalx, fthermaly
        end
        return fout, ftemp, feq, height, velx, vely, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py
    elseif device == "GPU"
        # Meso
        fout = CUDA.zeros(sysc.Lx, sysc.Ly, 9)
        ftemp = CUDA.zeros(sysc.Lx, sysc.Ly, 9)
        feq = CUDA.zeros(sysc.Lx, sysc.Ly, 9)
        
        # Macro
        height = CUDA.ones(sysc.Lx, sysc.Ly)
        velx = CUDA.zeros(sysc.Lx, sysc.Ly)
        vely = CUDA.zeros(sysc.Lx, sysc.Ly)
        pressure = zeros(sysc.Lx, sysc.Ly)

        # Forces
        Fx = CUDA.zeros(sysc.Lx, sysc.Ly)
        Fy = CUDA.zeros(sysc.Lx, sysc.Ly)
        slipx = CUDA.zeros(sysc.Lx, sysc.Ly)
        slipy = CUDA.zeros(sysc.Lx, sysc.Ly)
        h∇px = CUDA.zeros(sysc.Lx, sysc.Ly)
        h∇py = CUDA.zeros(sysc.Lx, sysc.Ly)
        if exotic
            fthermalx = CUDA.zeros(sysc.Lx, sysc.Ly)
            fthermaly = CUDA.zeros(sysc.Lx, sysc.Ly)
            return fout, ftemp, feq, height, velx, vely, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py, fthermalx, fthermaly
        end
        return fout, ftemp, feq, height, velx, vely, pressure, Fx, Fy, slipx, slipy, h∇px, h∇py
    end
end