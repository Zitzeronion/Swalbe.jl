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
  Expression: all(feq[:, :, 1] .≈ 1 - (2 / 3) * 0.01)
```

# References

- [Salmon](https://www.ingentaconnect.com/contentone/jmr/jmr/1999/00000057/00000003/art00005#)
- [Dellar](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.036309)
- [Peng et al.](https://onlinelibrary.wiley.com/doi/full/10.1002/fld.4726)

"""
function equilibrium!(feq, height, velocityx, velocityy, vsquare, gravity)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(feq) 
    # Some constants, gravity and weights
    g0 = 1.5 * gravity
    w1 = 1/9
    w5 = 1/36

    vsquare .= velocityx .* velocityx .+ velocityy .* velocityy 

    # Zeroth dist
    f0 .= height .* (1 .- 5/6 .* gravity .* height .- 2/3 .* vsquare)
    # First
    f1 .= w1 .* height .* (g0 .* height .+ 3 .* velocityx .+ 4.5 .* velocityx.^2 .- 1.5 .* vsquare)
    # Second
    f2 .= w1 .* height .* (g0 .* height .+ 3 .* velocityy .+ 4.5 .* velocityy.^2 .- 1.5 .* vsquare)
    # Third
    f3 .= w1 .* height .* (g0 .* height .- 3 .* velocityx .+ 4.5 .* velocityx.^2 .- 1.5 .* vsquare)
    # Forth
    f4 .= w1 .* height .* (g0 .* height .- 3 .* velocityy .+ 4.5 .* velocityy.^2 .- 1.5 .* vsquare)
    # Fifth
    f5 .= w5 .* height .* (g0 .* height .+ 3 .* (velocityx .+ velocityy) .+ 
                         4.5 .* (velocityx .+ velocityy).^2 .- 1.5 .* vsquare)
    # Sixth
    f6 .= w5 .* height .* (g0 .* height .+ 3 .* (velocityy .- velocityx) .+ 
                         4.5 .* (velocityy .- velocityx).^2 .- 1.5 .* vsquare)
    # Seventh
    f7 .= w5 .* height .* (g0 .* height .- 3 .* (velocityx .+ velocityy) .+
                         4.5 .* (velocityx .+ velocityy).^2 .- 1.5 .* vsquare)
    # Eigth
    f8 .= w5 .* height .* (g0 .* height .+ 3 .* (velocityx .- velocityy) .+ 
                         4.5 .* (velocityx .- velocityy).^2 .- 1.5 .* vsquare)
    return nothing
end

"""
    Equilibrium!(state, sys; g=sys.param.g)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann method.

# Arguments

- `state :: LBM_state_2D`: State data structure
- `sys :: SysConst`: System constants, contains the value for gravity
- `g <: Number`: Gravity, default is `sys.param.g`

# Examples
```jldoctest
julia> using Swalbe, Test

julia> sys = Swalbe.SysConst{Float64}(Lx=5, Ly=5, param=Swalbe.Taumucs(g=0.1)); state = Swalbe.Sys(sys, "CPU");

julia> state.velx .= fill(0.1,5,5);

julia> Swalbe.equilibrium!(state, sys) # Supply dummy u^2 as well.

julia> state.feq[:,:,1]
5×5 Matrix{Float64}:
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91
 0.91  0.91  0.91  0.91  0.91

julia> Swalbe.equilibrium!(state, sys, g=0.0) # Without gravity

julia> @test all(state.feq[:,:,1] .≈ 1 - 2/3 * 0.01)
Test Passed
  Expression: all(state.feq[:, :, 1] .≈ 1 - (2 / 3) * 0.01)
```

"""
function equilibrium!(state::LBM_state_2D, sys::SysConst; g=sys.param.g)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(state.feq) 
    # Some constants, gravity and weights
    g0 = 1.5 * g
    w1 = 1/9
    w5 = 1/36

    state.vsq .= state.velx .* state.velx .+ state.vely .* state.vely 

    # Zeroth dist
    f0 .= state.height .* (1 .- 5/6 .* g .* state.height .- 2/3 .* state.vsq)
    # First
    f1 .= w1 .* state.height .* (g0 .* state.height .+ 3 .* state.velx .+ 4.5 .* state.velx.^2 .- 1.5 .* state.vsq)
    # Second
    f2 .= w1 .* state.height .* (g0 .* state.height .+ 3 .* state.vely .+ 4.5 .* state.vely.^2 .- 1.5 .* state.vsq)
    # Third
    f3 .= w1 .* state.height .* (g0 .* state.height .- 3 .* state.velx .+ 4.5 .* state.velx.^2 .- 1.5 .* state.vsq)
    # Forth
    f4 .= w1 .* state.height .* (g0 .* state.height .- 3 .* state.vely .+ 4.5 .* state.vely.^2 .- 1.5 .* state.vsq)
    # Fifth
    f5 .= w5 .* state.height .* (g0 .* state.height .+ 3 .* (state.velx .+ state.vely) .+ 
                         4.5 .* (state.velx .+ state.vely).^2 .- 1.5 .* state.vsq)
    # Sixth
    f6 .= w5 .* state.height .* (g0 .* state.height .+ 3 .* (state.vely .- state.velx) .+ 
                         4.5 .* (state.vely .- state.velx).^2 .- 1.5 .* state.vsq)
    # Seventh
    f7 .= w5 .* state.height .* (g0 .* state.height .- 3 .* (state.velx .+ state.vely) .+
                         4.5 .* (state.velx .+ state.vely).^2 .- 1.5 .* state.vsq)
    # Eigth
    f8 .= w5 .* state.height .* (g0 .* state.height .+ 3 .* (state.velx .- state.vely) .+ 
                         4.5 .* (state.velx .- state.vely).^2 .- 1.5 .* state.vsq)
    return nothing
end

"""
    Equilibrium!(state, sys; g=sys.param.g)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann method.

# Arguments

- `state :: Expanded_2D`: State data structure
- `sys :: SysConst`: System constants, contains the value for gravity
- `g <: Number`: Gravity, default is `sys.param.g`

"""
function equilibrium!(state::Expanded_2D, sys::SysConst; g=sys.param.g)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(state.basestate.feq) 
    # Some constants, gravity and weights
    g0 = 1.5 * g
    w1 = 1/9
    w5 = 1/36

    state.basestate.vsq .= state.basestate.velx .* state.basestate.velx .+ state.basestate.vely .* state.basestate.vely 

    # Zeroth dist
    f0 .= state.basestate.height .* (1 .- 5/6 .* g .* state.basestate.height .- 2/3 .* state.basestate.vsq)
    # First
    f1 .= w1 .* state.basestate.height .* (g0 .* state.basestate.height .+ 3 .* state.basestate.velx .+ 4.5 .* state.basestate.velx.^2 .- 1.5 .* state.basestate.vsq)
    # Second
    f2 .= w1 .* state.basestate.height .* (g0 .* state.basestate.height .+ 3 .* state.basestate.vely .+ 4.5 .* state.basestate.vely.^2 .- 1.5 .* state.basestate.vsq)
    # Third
    f3 .= w1 .* state.basestate.height .* (g0 .* state.basestate.height .- 3 .* state.basestate.velx .+ 4.5 .* state.basestate.velx.^2 .- 1.5 .* state.basestate.vsq)
    # Forth
    f4 .= w1 .* state.basestate.height .* (g0 .* state.basestate.height .- 3 .* state.basestate.vely .+ 4.5 .* state.basestate.vely.^2 .- 1.5 .* state.basestate.vsq)
    # Fifth
    f5 .= w5 .* state.basestate.height .* (g0 .* state.basestate.height .+ 3 .* (state.basestate.velx .+ state.basestate.vely) .+ 
                         4.5 .* (state.basestate.velx .+ state.basestate.vely).^2 .- 1.5 .* state.basestate.vsq)
    # Sixth
    f6 .= w5 .* state.basestate.height .* (g0 .* state.basestate.height .+ 3 .* (state.basestate.vely .- state.basestate.velx) .+ 
                         4.5 .* (state.basestate.vely .- state.basestate.velx).^2 .- 1.5 .* state.basestate.vsq)
    # Seventh
    f7 .= w5 .* state.basestate.height .* (g0 .* state.basestate.height .- 3 .* (state.basestate.velx .+ state.basestate.vely) .+
                         4.5 .* (state.basestate.velx .+ state.basestate.vely).^2 .- 1.5 .* state.basestate.vsq)
    # Eight
    f8 .= w5 .* state.basestate.height .* (g0 .* state.basestate.height .+ 3 .* (state.basestate.velx .- state.basestate.vely) .+ 
                         4.5 .* (state.basestate.velx .- state.basestate.vely).^2 .- 1.5 .* state.basestate.vsq)
    return nothing
end

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
  Expression: all(feq[:, 1] .≈ (1 - 0.1 / 2) - 0.01)
```

"""
function equilibrium!(feq, height, velocity, gravity)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(feq) 
    
    # Zeroth dist
    f0 .= height .* (1 .- 0.5 .* gravity .* height .- velocity.^2)
    # First
    f1 .= height .* (0.25 .* gravity .* height .+ 0.5 .* velocity .+ 0.5 .* velocity.^2)
    # Second
    f2 .= height .* (0.25 .* gravity .* height .- 0.5 .* velocity .+ 0.5 .* velocity.^2)
    
    return nothing
end

"""
Equilibrium!(state, sys; g=sys.param.g)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann method.

# Arguments

- `state :: State_1D`: State data structure
- `sys :: SysConst`: System constants, contains the value for gravity
- `g <: Number`: Gravity, default is `sys.param.g`

"""
function equilibrium!(state::State_1D, sys::SysConst_1D; g=sys.param.g)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(state.feq) 
    
    # Zeroth dist
    f0 .= state.height .* (1 .- 0.5 .* g .* state.height .- state.vel.^2)
    # First
    f1 .= state.height .* (0.25 .* g .* state.height .+ 0.5 .* state.vel .+ 0.5 .* state.vel.^2)
    # Second
    f2 .= state.height .* (0.25 .* g .* state.height .- 0.5 .* state.vel .+ 0.5 .* state.vel.^2)
    
    return nothing
end

"""
Equilibrium!(state, sys; g=sys.param.g)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann method.

# Arguments

- `state :: Expanded_1D`: State data structure
- `sys :: SysConst`: System constants, contains the value for gravity
- `g <: Number`: Gravity, default is `sys.param.g`

"""
function equilibrium!(state::Expanded_1D, sys::SysConst_1D; g=sys.param.g)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(state.basestate.feq) 
    
    # Zeroth dist
    f0 .= state.basestate.height .* (1 .- 0.5 .* g .* state.basestate.height .- state.basestate.vel.^2)
    # First
    f1 .= state.basestate.height .* (0.25 .* g .* state.basestate.height .+ 0.5 .* state.basestate.vel .+ 0.5 .* state.basestate.vel.^2)
    # Second
    f2 .= state.basestate.height .* (0.25 .* g .* state.basestate.height .- 0.5 .* state.basestate.vel .+ 0.5 .* state.basestate.vel.^2)
    return nothing
end