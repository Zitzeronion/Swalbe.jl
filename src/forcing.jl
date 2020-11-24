"""
    slippage!(slipx, slipy, height, velx, vely, δ, μ)

Fluid substrate interaction that effectively mimics a velocity boundary condition at ``h=0``.

# Arguments

- `slipx::Array{Number,2}`: The x-component of the force due the velocity boundary condition
- `slipy::Array{Number,2}`: The y-component of the force due the velocity boundary condition
- `height::Array{Number,2}`: Height field `` h(\\mathbf{x},t)``
- `velx::Array{Number,2}`: x-component of the macroscopic velocity vector
- `vely::Array{Number,2}`: y-component of the macroscopic velocity vector
- `δ<:Number`: Extrapolation length into the substrate where the **no-slip** is met
- `μ<:Number`: Kinematic viscosity of the simulation, dependent on the value of **τ**

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



