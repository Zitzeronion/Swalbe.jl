"""
    Equilibrium!(feq, height,velocityx, velocityy, vsquare, gravity)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann method.

# Arguments

- `feq :: Array{<:Number,3}`: Equilibrium distribution function, to be calculated
- `height :: Array{<:Number,2}`: The height field `` h(\\mathbf{x},t)``
- `velocityx :: Array{<:Number,2}`: x-component of the macroscopic velocity
- `velocityy :: Array{<:Number,2}`: y-component of the macroscopic velocity
- `vsquare :: Array{<:Number,2}`: Dummy array that is preallocated to be filled with the square of the velocity vector
- `gravity <: Number`: Strength of the gravitational acceleration in lattice units

# Mathematics

The detailed motivation and derivation can be found in the article of Salmon.
Similar to the standard Navier-Stokes approximating lattice Boltzmann methods the equilibrium distribution `feq`
is an expansion to second order in velocity `u`.
If you want it is an slow mode part of the shallow water theory and thus the equilibrium is given as

`` f_i^{\\text{eq}} = h \\bigg(1 - \\frac{5}{6}g h - \\frac{2}{3}u^2\\bigg),\\quad i=0 \\newline
   f_i^{\\text{eq}} = w_i h \\bigg(g h + 3 \\mathbf{c}_i\\cdot\\mathbf{u} + \\frac{9}{2}(\\mathbf{c}_i\\cdot\\mathbf{u})^2) + \\frac{3}{2}u^2\\bigg),\\quad else ``

where ``g`` is the gravitational acceleration (in lattice units) and ``w_i, \\mathbf{c}_i`` are the weights and lattice velocities.

It has been shown that it is possible to get rid of the gravity driven term in the equilibrium distribution, thus

`` f_i^{\\text{eq}} = h \\bigg(1 - \\frac{2}{3}u^2\\bigg),\\quad i=0 \\newline
f_i^{\\text{eq}} = w_i h \\bigg(3 \\mathbf{c}_i\\cdot\\mathbf{u} + \\frac{9}{2}(\\mathbf{c}_i\\cdot\\mathbf{u})^2) + \\frac{3}{2}u^2\\bigg),\\quad else ``

then of course the topography gradient has to be included as a force term.

# Examples
```jldoctest
julia> using Swalbe, Test

julia> feq = zeros(5,5,9); ρ = ones(5,5); ux = fill(0.1,5,5); uy = zeros(5,5);

julia> Swalbe.equilibrium!(feq, ρ, ux, uy, zeros(5,5), 0.1) # Supply dummy u^2 as well.

julia> feq[:,:,1]
5×5 Matrix{Float64}:
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91

julia> Swalbe.equilibrium!(feq, ρ, ux, uy, zeros(5,5), 0.0) # Supply dummy u^2 as well.

julia> @test all(feq[:,:,1] .≈ 1 - 2/3 * 0.01)
Test Passed
```

# References

- [Salmon](https://www.ingentaconnect.com/contentone/jmr/jmr/1999/00000057/00000003/art00005#)
- [Dellar](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.036309)
- [Peng et al.](https://onlinelibrary.wiley.com/doi/full/10.1002/fld.4726)

"""
function equilibrium!(feq, height, velocityx, velocityy, vsquare, gravity)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(feq) 
    # Gravity and weights
    g0 = 1.5 * gravity
    w1 = 1/9
    w5 = 1/36
    # Expansion coefficients velocity space
    b = 3
    c = 4.5
    d = 1.5

    vsquare .= @. velocityx * velocityx + velocityy * velocityy 

    # Zeroth dist
    f0 .= @. height * (1 - 5/6 * gravity * height - 2/3 * vsquare)
    # First
    f1 .= @. w1 * height * (g0 * height + b * velocityx + c * velocityx^2 - d * vsquare)
    # Second
    f2 .= @. w1 * height * (g0 * height + b * velocityy + c * velocityy^2 - d * vsquare)
    # Third
    f3 .= @. w1 * height * (g0 * height - b * velocityx + c * velocityx^2 - d * vsquare)
    # Forth
    f4 .= @. w1 * height * (g0 * height - b * velocityy + c * velocityy^2 - d * vsquare)
    # Fifth
    f5 .= @. w5 * height * (g0 * height + b * (velocityx + velocityy) + c * (velocityx + velocityy)^2 - d * vsquare)
    # Sixth
    f6 .= @. w5 * height * (g0 * height + b * (velocityy - velocityx) + c * (velocityy - velocityx)^2 - d * vsquare)
    # Seventh
    f7 .= @. w5 * height * (g0 * height - b * (velocityx + velocityy) + c * (velocityx + velocityy)^2 - d * vsquare)
    # Eigth
    f8 .= @. w5 * height * (g0 * height + b * (velocityx - velocityy) + c * (velocityx - velocityy)^2 - d * vsquare)
    return nothing
end

equilibrium!(state::LBM_state_2D, sys::SysConst) = equilibrium!(state.feq, state.height, state.velx, state.vely, state.vsq, sys.param.g)

equilibrium!(state::Expanded_2D, sys::SysConst) = equilibrium!(state.basestate.feq, state.basestate.height, state.basestate.velx, state.basestate.vely, state.basestate.vsq, sys.param.g)


"""
Equilibrium!(feq, height, velocity, gravity)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann method.

# Arguments

- `feq :: Array{<:Number,3}`: Equilibrium distribution function, to be calculated
- `height :: Array{<:Number,2}`: The height field `` h(\\mathbf{x},t)``
- `velocity :: Array{<:Number,2}`: x-component of the macroscopic velocity
- `gravity <: Number`: Strength of the gravitational acceleration in lattice units

# Examples
```jldoctest
julia> using Swalbe, Test

julia> feq = zeros(10,3); ρ = ones(10); u = fill(0.1,10);

julia> Swalbe.equilibrium!(feq, ρ, u, 0.1) # Supply dummy u^2 as well.

julia> feq[:,1]
10-element Vector{Float64}:
 0.94
 0.94
 0.94
 0.94
 0.94
 0.94
 0.94
 0.94
 0.94
 0.94

julia> @test all(feq[:,1] .≈ 1 - 0.1/2 - 0.01)
Test Passed
```

"""
function equilibrium!(feq, height, velocity, gravity)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(feq) 
    a = 0.25
    b = 0.5

    # Zeroth dist
    f0 .= @. height * (1 - b * gravity * height - velocity^2)
    # First
    f1 .= @. height * (a * gravity * height + b * velocity + b * velocity^2)
    # Second
    f2 .= @. height * (a * gravity * height - b * velocity + b * velocity^2)
    
    return nothing
end

equilibrium!(state::State_1D, sys::Consts_1D; g=sys.param.g) = equilibrium!(state.feq, state.height, state.vel, sys.param.g)

equilibrium!(state::Expanded_1D, sys::Consts_1D) = equilibrium!(state.basestate.feq, state.basestate.height, state.basestate.vel, sys.param.g)
