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

julia> Swalbe.equilibrium!(feq, ρ, ux, uy, zeros(5,5)) # Use the dispatch without g.

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

# Dispatch without gravity input
function equilibrium!(feq, height, velocityx, velocityy, vsquare)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(feq) 
    # Some constants, gravity and weights
    w1 = 1/9
    w5 = 1/36

    vsquare .= velocityx .* velocityx .+ velocityy .* velocityy 

    # Zeroth dist
    f0 .= height .* (1 .- 2/3 .* vsquare)
    # First
    f1 .= w1 .* height .* (3 .* velocityx .+ 4.5 .* velocityx.^2 .- 1.5 .* vsquare)
    # Second
    f2 .= w1 .* height .* (3 .* velocityy .+ 4.5 .* velocityy.^2 .- 1.5 .* vsquare)
    # Third
    f3 .= w1 .* height .* (-3 .* velocityx .+ 4.5 .* velocityx.^2 .- 1.5 .* vsquare)
    # Forth
    f4 .= w1 .* height .* (-3 .* velocityy .+ 4.5 .* velocityy.^2 .- 1.5 .* vsquare)
    # Fifth
    f5 .= w5 .* height .* (3 .* (velocityx .+ velocityy) .+ 
                         4.5 .* (velocityx .+ velocityy).^2 .- 1.5 .* vsquare)
    # Sixth
    f6 .= w5 .* height .* (3 .* (velocityy .- velocityx) .+ 
                         4.5 .* (velocityy .- velocityx).^2 .- 1.5 .* vsquare)
    # Seventh
    f7 .= w5 .* height .* (-3 .* (velocityx .+ velocityy) .+
                         4.5 .* (velocityx .+ velocityy).^2 .- 1.5 .* vsquare)
    # Eigth
    f8 .= w5 .* height .* (3 .* (velocityx .- velocityy) .+ 
                         4.5 .* (velocityx .- velocityy).^2 .- 1.5 .* vsquare)
    return nothing
end

function equilibrium!(state::State, sys::SysConst)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(state.feq) 
    # Some constants, gravity and weights
    g0 = 1.5 * sys.g
    w1 = 1/9
    w5 = 1/36

    state.vsq .= state.velx .* state.velx .+ state.vely .* state.vely 

    # Zeroth dist
    f0 .= state.height .* (1 .- 5/6 .* sys.g .* state.height .- 2/3 .* state.vsq)
    # First
    f1 .= w1 .* state.height .* (g0 .* state.height .+ 3 .* state.velx .+ 4.5 .* state.velx .* state.velx .- 1.5 .* state.vsq)
    # Second
    f2 .= w1 .* state.height .* (g0 .* state.height .+ 3 .* state.vely .+ 4.5 .* state.vely .* state.vely .- 1.5 .* state.vsq)
    # Third
    f3 .= w1 .* state.height .* (g0 .* state.height .- 3 .* state.velx .+ 4.5 .* state.velx .* state.velx .- 1.5 .* state.vsq)
    # Forth
    f4 .= w1 .* state.height .* (g0 .* state.height .- 3 .* state.vely .+ 4.5 .* state.vely .* state.vely .- 1.5 .* state.vsq)
    # Fifth
    f5 .= w5 .* state.height .* (g0 .* state.height .+ 3 .* (state.velx .+ state.vely) .+ 
                         4.5 .* (state.velx .+ state.vely) .* (state.velx .+ state.vely) .- 1.5 .* state.vsq)
    # Sixth
    f6 .= w5 .* state.height .* (g0 .* state.height .+ 3 .* (state.vely .- state.velx) .+ 
                         4.5 .* (state.vely .- state.velx) .* (state.vely .- state.velx) .- 1.5 .* state.vsq)
    # Seventh
    f7 .= w5 .* state.height .* (g0 .* state.height .- 3 .* (state.velx .+ state.vely) .+
                         4.5 .* (state.velx .+ state.vely) .* (state.velx .+ state.vely) .- 1.5 .* state.vsq)
    # Eigth
    f8 .= w5 .* state.height .* (g0 .* state.height .+ 3 .* (state.velx .- state.vely) .+ 
                         4.5 .* (state.velx .- state.vely) .* (state.velx .- state.vely) .- 1.5 .* state.vsq)
    return nothing
end
# Dispatch without gravity and state struct
function equilibrium!(state::State)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(state.feq) 
    # Some constants, gravity and weights
    w1 = 1/9
    w5 = 1/36

    state.vsq .= state.velx .* state.velx .+ state.vely .* state.vely 

    # Zeroth dist
    f0 .= state.height .* (1 .- 2/3 .* state.vsq)
    # First
    f1 .= w1 .* state.height .* (3 .* state.velx .+ 4.5 .* state.velx .* state.velx .- 1.5 .* state.vsq)
    # Second
    f2 .= w1 .* state.height .* (3 .* state.vely .+ 4.5 .* state.vely .* state.vely .- 1.5 .* state.vsq)
    # Third
    f3 .= w1 .* state.height .* (-3 .* state.velx .+ 4.5 .* state.velx .* state.velx .- 1.5 .* state.vsq)
    # Forth
    f4 .= w1 .* state.height .* (-3 .* state.vely .+ 4.5 .* state.vely .* state.vely .- 1.5 .* state.vsq)
    # Fifth
    f5 .= w5 .* state.height .* (3 .* (state.velx .+ state.vely) .+ 
                         4.5 .* (state.velx .+ state.vely) .* (state.velx .+ state.vely) .- 1.5 .* state.vsq)
    # Sixth
    f6 .= w5 .* state.height .* (3 .* (state.vely .- state.velx) .+ 
                         4.5 .* (state.vely .- state.velx) .* (state.vely .- state.velx) .- 1.5 .* state.vsq)
    # Seventh
    f7 .= w5 .* state.height .* (-3 .* (state.velx .+ state.vely) .+
                         4.5 .* (state.velx .+ state.vely) .* (state.velx .+ state.vely) .- 1.5 .* state.vsq)
    # Eight
    f8 .= w5 .* state.height .* (3 .* (state.velx .- state.vely) .+ 
                         4.5 .* (state.velx .- state.vely) .* (state.velx .- state.vely) .- 1.5 .* state.vsq)
    return nothing
end
#active thin film no gravity
function equilibrium!(state::StateActive)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = viewdists(state.feq) 
    # Some constants, gravity and weights
    w1 = 1/9
    w5 = 1/36

    state.vsq .= state.velx .* state.velx .+ state.vely .* state.vely 

    # Zeroth dist
    f0 .= state.height .* (1 .- 2/3 .* state.vsq)
    # First
    f1 .= w1 .* state.height .* (3 .* state.velx .+ 4.5 .* state.velx .* state.velx .- 1.5 .* state.vsq)
    # Second
    f2 .= w1 .* state.height .* (3 .* state.vely .+ 4.5 .* state.vely .* state.vely .- 1.5 .* state.vsq)
    # Third
    f3 .= w1 .* state.height .* (-3 .* state.velx .+ 4.5 .* state.velx .* state.velx .- 1.5 .* state.vsq)
    # Forth
    f4 .= w1 .* state.height .* (-3 .* state.vely .+ 4.5 .* state.vely .* state.vely .- 1.5 .* state.vsq)
    # Fifth
    f5 .= w5 .* state.height .* (3 .* (state.velx .+ state.vely) .+ 
                         4.5 .* (state.velx .+ state.vely) .* (state.velx .+ state.vely) .- 1.5 .* state.vsq)
    # Sixth
    f6 .= w5 .* state.height .* (3 .* (state.vely .- state.velx) .+ 
                         4.5 .* (state.vely .- state.velx) .* (state.vely .- state.velx) .- 1.5 .* state.vsq)
    # Seventh
    f7 .= w5 .* state.height .* (-3 .* (state.velx .+ state.vely) .+
                         4.5 .* (state.velx .+ state.vely) .* (state.velx .+ state.vely) .- 1.5 .* state.vsq)
    # Eight
    f8 .= w5 .* state.height .* (3 .* (state.velx .- state.vely) .+ 
                         4.5 .* (state.velx .- state.vely) .* (state.velx .- state.vely) .- 1.5 .* state.vsq)
    return nothing
end

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

function equilibrium!(feq, height, velocity)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(feq) 
    
    # Zeroth dist
    f0 .= height .* (1 .- velocity.^2)
    # First
    f1 .= height .* (0.5 .* velocity .+ 0.5 .* velocity.^2)
    # Second
    f2 .= height .* (-0.5 .* velocity .+ 0.5 .* velocity.^2)
    return nothing
end

function equilibrium!(state::State_1D, sys::SysConst_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(state.feq) 
    
    # Zeroth dist
    f0 .= state.height .* (1 .- 0.5 .* sys.g .* state.height .- state.vel .* state.vel)
    # First
    f1 .= state.height .* (0.25 .* sys.g .* state.height .+ 0.5 .* state.vel .+ 0.5 .* state.vel .* state.vel)
    # Second
    f2 .= state.height .* (0.25 .* sys.g .* state.height .- 0.5 .* state.vel .+ 0.5 .* state.vel .* state.vel)
    
    return nothing
end

function equilibrium!(state::State_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(state.feq) 
    
    # Zeroth dist
    f0 .= state.height .* (1 .- state.vel .* state.vel)
    # First
    f1 .= state.height .* (0.5 .* state.vel .+ 0.5 .* state.vel .* state.vel)
    # Second
    f2 .= state.height .* (-0.5 .* state.vel .+ 0.5 .* state.vel .* state.vel)
    return nothing
end


function equilibrium_new!(state::State_1D, sys::SysConst_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(state.feq) 
    
    # Zeroth dist
    f0 .= state.height 
    # First
    f1 .=(-0.5 .*( state.h∇p ./ state.slip ))
    # Second
    f2 .= (0.5 .* (state.h∇p ./ state.slip ))
    return nothing
end
#thin film no gravity
function equilibrium!(state::Active_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    f0, f1, f2 = viewdists_1D(state.feq) 
    
    # Zeroth dist
    f0 .= state.height .* (1 .- state.vel .* state.vel)
    # First
    f1 .= state.height .* (0.5 .* state.vel .+ 0.5 .* state.vel .* state.vel)
    # Second
    f2 .= state.height .* (-0.5 .* state.vel .+ 0.5 .* state.vel .* state.vel)
    return nothing
end


"""
    rho_equilibrium_quadratic!(state)
Calculates the equilibrium for ``\\rho`` exactly as for h.
"""
function rho_equilibrium_quadratic!(state::Active_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    g0, g1, g2 = viewdists_1D(state.geq) 
    
    # Zeroth dist
    g0 .= state.rho .* (1 .- state.rho_vel .* state.rho_vel)
    # First
    g1 .= state.rho .* (0.5 .* state.rho_vel .+ 0.5 .* state.rho_vel .* state.rho_vel)
    # Second
    g2 .= state.rho .* (-0.5 .* state.rho_vel .+ 0.5 .* state.rho_vel .* state.rho_vel)
    return nothing
end

"""
    rho_A_equilibrium_quadratic!(state)
Calculates the equilibrium for ``\\rho`` exactly as for h.
"""
function rho_A_equilibrium_quadratic!(state::Active_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    g0, g1, g2 = viewdists_1D(state.heq) 
    
    # Zeroth dist
    g0 .= state.rho_A .* (1 .- state.rho_A_vel .* state.rho_A_vel)
    # First
    g1 .= state.rho_A .* (0.5 .* state.rho_A_vel .+ 0.5 .* state.rho_A_vel .* state.rho_A_vel)
    # Second
    g2 .= state.rho_A .* (-0.5 .* state.rho_A_vel .+ 0.5 .* state.rho_A_vel .* state.rho_A_vel)
    return nothing
end



"""
    rho_B_equilibrium_quadratic!(state)
Calculates the equilibrium for ``\\rho`` exactly as for h.
"""
function rho_B_equilibrium_quadratic!(state::Active_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    g0, g1, g2 = viewdists_1D(state.beq) 
    
    # Zeroth dist
    g0 .= state.rho_B .* (1 .- state.rho_B_vel .* state.rho_B_vel)
    # First
    g1 .= state.rho_B .* (0.5 .* state.rho_B_vel .+ 0.5 .* state.rho_B_vel .* state.rho_B_vel)
    # Second
    g2 .= state.rho_B .* (-0.5 .* state.rho_B_vel .+ 0.5 .* state.rho_B_vel .* state.rho_B_vel)
    return nothing
end



"""
	function rho_equilibrium_linear!(state::Active_1D, sys::SysConstActive_1D)

Linear equilibrium, untested and unused.
"""
function rho_equilibrium_linear!(state::Active_1D, sys::SysConstActive_1D)
    # Views help to circumvent having a loop, which sucks on the GPU
    g0, g1, g2 = viewdists_1D(state.geq) 

    # calculate derivatives as finite differences
    ∇f!(state.grad_rho, state.rho, state.dgrad)
    ∇f!(state.grad_h, state.height, state.dgrad)
    # Zeroth dist
    g0 .= state.rho 
    # First
    # Another advective term has to be added here if desired
    g1 .=0.5.* (state.rho   .* ( sys.D .* state.grad_h ./ (state.height.-sys.R)) .+ state.vel - sys.D .* state.grad_rho)
    # Second
    g2 .=0.5.* (state.rho   .* ( .- sys.D .*state.grad_h ./ (state.height.- sys.R)) .- state.vel.+ sys.D .* state.grad_rho)
    return nothing
end
