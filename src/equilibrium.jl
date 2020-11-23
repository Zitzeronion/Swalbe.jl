"""
    Equilibrium!(feq, height,velocityx, velocityy, vsquare, gravity)

Calculation of the equilibrium distribution `feq` for the shallow water lattice Boltzmann.

# Mathematics

The detailed motivation and derivation can be found in the article of Salmon.
Similar to the standard Navier-Stokes approximating lattice Boltzmann methods the equilibrium distribution `feq`
is an expansion to second order in velocity `u`.
If you want it is an slow mode part of the shallow water theory and thus the equilibrium is given as

`` f_i^{\\text{eq}} = h \\bigg(1 - \\frac{5}{6}g h - \\frac{2}{3}u^2\\bigg),\\quad i=0 \\newline
   f_i^{\\text{eq}} = w_i h \\bigg(g h + 3 \\mathbf{c}_i\\cdot\\mathbf{u} + \\frac{9}{2}(\\mathbf{c}_i\\cdot\\mathbf{u})^2) + \\frac{3}{2}u^2,\\quad else ``

where ``g`` is the strength of the gravitational acceleration and ``w_i, \\mathbf{c}_i`` are the weights and lattice velocities. 
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