"""
    slippage!(slipx, slipy, height, velx, vely, δ, μ)

Fluid substrate interaction that effectively mimics a velocity boundary condition at ``h=0``.

# Arguments

- `slipx :: Array{<:Number,2}`: The x-component of the force due the velocity boundary condition
- `slipy :: Array{<:Number,2}`: The y-component of the force due the velocity boundary condition
- `height::Array{<:Number,2}`: Height field `` h(\\mathbf{x},t)``
- `velx::Array{<:Number,2}`: x-component of the macroscopic velocity vector
- `vely::Array{<:Number,2}`: y-component of the macroscopic velocity vector
- `δ <: Number`: Extrapolation length into the substrate where the **no-slip** is met
- `μ <: Number`: Kinematic viscosity of the simulation, dependent on the value of **τ**

# Mathematics

With the velocity boundary condition at the fluid substrate interface we build the second main descriptor between our model and the thin film equation.
One well studied assumption is that the fluid velocity vanishes at `` h(\\mathbf{x}) = 0 `` which is called **no slip** condition.

TBD

# Examples

TBD

# References

- [Zitz, Scagliarini and Harting](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)
- [Münch, Wagner und Witelski](https://link.springer.com/article/10.1007/s10665-005-9020-3)
- [Oron, Davis and Bankoff](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.69.931)
"""
function slippage!(slipx, slipy, height, velx, vely, δ, μ)
    slipx .= (6μ .* height .* velx) ./ (2 .* height.^2 .+ 6δ .* height .+ 3δ^2 )
    slipy .= (6μ .* height .* vely) ./ (2 .* height.^2 .+ 6δ .* height .+ 3δ^2 )
    return nothing
end

"""
    thermal!(fx, fy, height, kᵦT, μ, δ)

Computations of force due to thermal fluctuations.

# Arguments

- `fx :: Array{<:Number,2}`: x-component of the force due to the fluctuations
- `fy :: Array{<:Number,2}`: y-component of the force due to the fluctuations
- `height :: Array{<:Number,2}`: Height field ``h(\\mathbf{x},t) ``
- `kᵦT <: Number`: Strenght of thermal fluctuations in lattice units
- `μ <: Number`: The kinetic viscosity
- `δ <: Number`: Slip length, needed to normalize

# Mathematics

TBD

# Examples
```jldoctest
julia> using Swalbe, Statistics, Test

julia> x = ones(50,50); y = ones(50,50); h = ones(50,50);

julia> Swalbe.thermal!(x, y, h, 0.1, 1/6, 1)

julia> @test mean(x) ≈ 0.0 atol=1e-2
Test Passed

julia> @test mean(y) ≈ 0.0 atol=1e-2
Test Passed

julia> @test var(x) ≈ 0.0 atol=1e-2
Test Passed

```

# References

- [Grün, Mecke and Rauscher](https://link.springer.com/article/10.1007/s10955-006-9028-8)
- [Mecke, Rauscher](https://iopscience.iop.org/article/10.1088/0953-8984/17/45/042/meta)
- [Davidovitch, Moro and Stone](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.244505)

""" 
function thermal!(fluc_x::CuArray, fluc_y::CuArray, height, kᵦT , μ, δ)
    len, wid = size(height)
    fluc_x .= CUDA.randn(len, wid)
    fluc_y .= CUDA.randn(len, wid)
    fluc_x .*= sqrt.(2 .* kᵦT .* μ .* 6 .* height ./
                    (2 .* height.^2 .+
                     6 .* height .* δ .+
                     3 .* δ^2))
    fluc_y .*= sqrt.(2 .* kᵦT .* μ .* 6 .* height ./
                    (2 .* height.^2 .+
                     6 .* height .* δ .+
                     3 .* δ^2))
    return nothing
end
# Apparently one can not use map(rand, array) with CuArrays
function thermal!(fluc_x::Array, fluc_y::Array, height, kᵦT , μ, δ)
    len, wid = size(height)
    fluc_x .= randn(len, wid)
    fluc_y .= randn(len, wid)
    fluc_x .*= sqrt.(2 .* kᵦT .* μ .* 6 .* height ./
                    (2 .* height.^2 .+
                     6 .* height .* δ .+
                     3 .* δ^2))
    fluc_y .*= sqrt.(2 .* kᵦT .* μ .* 6 .* height ./
                    (2 .* height.^2 .+
                     6 .* height .* δ .+
                     3 .* δ^2))
    return nothing
end
