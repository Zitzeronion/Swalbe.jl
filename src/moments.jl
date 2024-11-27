"""
    moments!(height, velx, vely, fout)

Computation of the hydrodynamic moments, `height` and `velocity`.

# Mathematics

The macroscopic quantities such as the height and the velocity are the moments of the distribution function, 

`` h(\\mathbf{x},t) = \\sum_{i=0}^8 f_i , ``

and 

`` \\mathbf{v}(\\mathbf{x},t) = \\frac{1}{h}\\sum_{i=0}^8 \\mathbf{c}_i f_i . ``

# Examples

```jldoctest
julia> using Swalbe, Test

julia> fout = zeros(5,5,9); fout[:,:,1] .= 1.0; fout[:,:,2] .= 0.1; # Dist with artifical velocity in x

julia> height = zeros(5,5); velx = zeros(5,5); vely = zeros(5,5);

julia> Swalbe.moments!(height,velx,vely,fout)

julia> @test all(height .== 1.1)
Test Passed

julia> @test all(velx .== 0.1/1.1)
Test Passed

julia> @test all(vely .== 0.0)
Test Passed
```

# References

- [Krüger](https://www.springer.com/gp/book/9783319446479)
- [Salmon](https://www.ingentaconnect.com/contentone/jmr/jmr/1999/00000057/00000003/art00005#)
- [Zitz, Scagliarini and Harting](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)
"""
function moments!(height, velx, vely, fout)
    # Get views of the populations
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(fout)
    # Compute the height
    sum!(height, fout)
    # and the velocities (as simple as possible)
    velx .= (f1 .- f3 .+ f5 .- f6 .- f7 .+ f8) ./ height
    vely .= (f2 .- f4 .+ f5 .+ f6 .- f7 .- f8) ./ height
    return nothing
end

function moments!(height::Vector, vel, fout)
    # Get views of the populations
    f0, f1, f2 = Swalbe.viewdists_1D(fout)
    # Compute the height
    sum!(height, fout)
    # and the velocities (as simple as possible)
    vel .= (f1 .- f2) ./ height
    return nothing
end

moments!(state::LBM_state_2D) = moments!(state.height, state.velx, state.vely, state.fout)
moments!(state::CuState) = moments!(state.height, state.velx, state.vely, state.fout)
moments!(state::LBM_state_1D) = moments!(state.height, state.vel, state.fout)

moments!(state::Expanded_2D) = moments!(
    state.basestate.height,
    state.basestate.velx,
    state.basestate.vely,
    state.basestate.fout,
)

moments!(state::Expanded_1D) = moments!(state.basestate.height, state.basestate.vel, state.basestate.fout)
moments!(state::Active_1D) = moments!(state.height, state.vel, state.fout)


"""
	function moments!(state::MultiLayer_2D, sys::SysConstMultiLayer)

Same as [moments!](@ref) but for multlayer calculating all layer height 
"""
function moments!(state::MultiLayer_2D, sys::SysConstMultiLayer)
    # Get views of the populations
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdistsMultiLayer(state.fout) 
    # Compute the height
    # I guess we are loosing some efficiency here by not doing all at once
    for i in 1:sys.layers
        height_view=view(state.height,:,:,i)
        sum!(height_view, state.fout[:,:,:,i])
    end
    # and the velocities (as simple as possible)
    state.velx .= (f1 .- f3 .+ f5 .- f6 .- f7 .+ f8) ./ state.height
    state.vely .= (f2 .- f4 .+ f5 .+ f6 .- f7 .- f8) ./ state.height
    return nothing
end


"""
    rho_moments!(state)

Calculates the moments for ``\\rho`` as
`` \\rho = \\sum g_i``
`` \\rho u = v - D\\frac{\\nabla \\rho}{\\rho}- DA_2 \\nabla \\rho + D \\nabla h \\nabla_h \\mathcal{F}``

That leads to solving
`` \\partial_t \\rho = \\nabla \\cdot\\left[ D(1 + A_2 \\rho) \\nabla  \\rho - \\rho u\\right]
``A_2`` is a second order diffusive term that standarized is choosen ``A_2=0``. Also we iclude a window integration to prevent numerical noise at step functions

`` \\rho(x) \\mapsto \\rho \\omega+ \\frac{1-\\omega}{2}(\\rho(x+\\Delta x) + \\rho(x-\\Delta x)``

where standarized `` \\omega =0``.

That is like the standard model for advection difusion as seen in [Tim Krueger Lattice Boltzman methods] but without the diffusive term in the equilibrium function.

See also: [`Swalbe.rho_equilibrium`](@ref)
"""

function rho_moments!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho, state.gout)
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho, state.rho, state.dgrad)
     ∇f!(state.grad_p, state.pressure, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
     state.D .= sys.D .* ( 1 .+ sys.A_2 .* state.rho)
     if sys.alpha != 0
        state.rho_adv_vel .= (
            state.∇γ ./ sys.μ .* state.V_gamma
            .- state.grad_p ./ sys.μ .* state.V_p
        )
       state.rho_vel .= state.rho_adv_vel  .- state.D .* state.grad_rho ./ state.rho  .- sys.D .* state.grad_h .* state.nabla_F
     else
        state.rho_vel .= state.vel .- state.D .* state.grad_rho ./ state.rho .+ sys.D .* state.grad_h ./ (state.height)
     end
     return nothing
 end

"""
    rho_A_moments!(state)

Calculates the moments for ``\\rho_A`` as

`` \\rho_A = \\sum h_i``

`` \\rho_A u = v - D_A\\frac{\\nabla \\rho_A}{\\rho_A}+ D_A \\nabla h \\nabla_h \\mathcal{F}``

That leads to solving

`` \\partial_t \\rho_A = \\nabla \\cdot\\left[ D_A \\nabla \\rho_A - \\rho_A u \\right]

We only have implemented ``\\alpha=0`` for ``\\rho_A`` i.e. we assume the product to be neutrally buoyant.

That is like the standard model for advection difusion as seen in [Tim Krueger Lattice Boltzman methods] but without the diffusive term in the equilibrium function.

See also: [`Swalbe.rho_A_equilibrium`](@ref)
"""

function rho_A_moments!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho_A, state.hout)
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho_A, state.rho_A, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
     state.rho_A_vel .= state.vel .- sys.D_Product .* state.grad_rho_A ./ state.rho_A .+ sys.D_Product .* state.grad_h ./ (state.height)
     return nothing
 end


"""
    rho_B_moments!(state::Active_1D, sys::SysConstActive_1D)

Calculates the moments for the secondary density component ``\\rho_B`` based on the sum of the particle distribution function components, ``b_i``:

`` \\rho_B = \\sum b_i``

The flux term ``\\rho_B u`` is defined to include advection and diffusive effects without buoyancy forces:

`` \\rho_B u = v - D_B \\frac{\\nabla \\rho_B}{\\rho_B} + D_B \\nabla h \\nabla_h \\mathcal{F}``

This formulation leads to an advection-diffusion equation:

`` \\partial_t \\rho_B = \\nabla \\cdot \\left[ D_B \\nabla \\rho_B - \\rho_B u \\right]``

where ``D_B`` is the diffusion coefficient specific to ``\\rho_B``. The model applies diffusion and advection similar to lattice Boltzmann methods (see [Tim Krueger, Lattice Boltzmann Methods]).

### Arguments
- `state::Active_1D`: The system's state, containing fields for densities, gradients, and velocities.
- `sys::SysConstActive_1D`: System constants for the diffusion and advection terms, specific to ``\\rho_B``.

### Returns
- `nothing`: The function modifies `state` in place, updating ``\\rho_B`` and its related velocity fields.

### See Also
- `rho_A_moments!`: A similar function for calculating moments for ``\\rho_A``.
"""

function rho_B_moments!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho_B, state.bout)
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho_B, state.rho_B, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
     state.rho_B_vel .= state.vel .- sys.D_B .* state.grad_rho_B ./ state.rho_B .+ sys.D_B .* state.grad_h ./ (state.height)
     return nothing
 end

 """
    moments!(state::StateMiscible_1D, sys::SysConstMiscible_1D)

Compute the hydrodynamic moments, including height, velocity, and diffusion effects for a miscible system in 1D.

# Description

This function calculates the macroscopic hydrodynamic quantities for a miscible fluid system in one dimension. Specifically, it computes:
1. The height (density) for each liquid component.
2. The velocity, considering contributions from advection and diffusion effects.

Additionally, it evaluates diffusion gradients and incorporates these into the velocity calculations using a miscible fluids model.

# Mathematics

For each liquid component \\( i \\), the height is computed as the zeroth moment of the distribution function:
```math
h_i(x, t) = \\sum_{j=0}^2 f_{ij},
```
where \\( f_{ij} \\) represents the distribution function for the \\( j \\)-th velocity direction of liquid \\( i \\).

The velocity is computed as:
```math
v(x, t) = \\frac{f_1 - f_2}{h},
```
where \\( f_1 \\) and \\( f_2 \\) are the populations corresponding to the positive and negative velocity directions, respectively.

Diffusion corrections are applied to the velocity based on the gradient of the height \\( h \\):
```math
v_i(x, t) -= D \\cdot \\frac{h_j \\cdot \\nabla h_i - h_i \\cdot \\nabla h_j}{(h_i + h_j) \\cdot h_i},
```
where \\( D \\) is the diffusion coefficient, and \\( \\nabla h \\) is the height gradient.


# Examples

```jldoctest
julia> using Swalbe, Test

julia> state = StateMiscible_1D(fout, height, vel, grad_h, dgrad); # Initialize state

julia> sys = SysConstMiscible_1D(D=0.01, liquids=2); # Define system constants

julia> fout[:,:,1] .= 1.0; fout[:,:,2] .= 0.5; # Set distribution function

julia> height = zeros(10, 2); vel = zeros(10, 2); grad_h = zeros(10, 2); dgrad = 1.0;

julia> Swalbe.moments!(state, sys)

julia> @test all(state.height[:,1] .== 1.5) # First liquid height
Test Passed

julia> @test all(state.vel[:,1] .== 0.01) # Corrected velocity for first liquid
Test Passed
```

# Arguments

- `state::StateMiscible_1D`: A structure holding the state of the miscible fluid system, including distribution functions, heights, velocities, and gradients.
- `sys::SysConstMiscible_1D`: A structure defining system-wide constants such as diffusion coefficient and number of liquid components.

# Returns

`nothing`. The results are stored in-place in the provided `state`.

# References

- [Krüger et al., *The Lattice Boltzmann Method*](https://www.springer.com/gp/book/9783319446479)
- [Salmon, *Hydrodynamics of Lattice Boltzmann*](https://www.ingentaconnect.com/contentone/jmr/jmr/1999/00000057/00000003/art00005#)
- [Zitz et al., *Multicomponent LBM Simulation*](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)
- [Baumgartner et al.](https://doi.org/10.1073/pnas.2120432119)
"""
function moments!(state::StateMiscible_1D, sys::SysConstMiscible_1D)
    # Get views of the populations
    f0, f1, f2 = Swalbe.viewdistsMiscible_1D(state.fout)
    # Compute the height
    for i in 1:sys.liquids
        height_view=view(state.height,:,i)
        sum!(height_view, state.fout[:,:,i])
    end
    # and the velocities (as simple as possible)
    state.vel .= (f1 .- f2) ./ state.height
    #diffusion
    ∇f_Miscible!(state.grad_h, state.height, state.dgrad)
    # state.vel[:,1] .-= sys.D .* ( (state.height[:,2] .* state.grad_h[:,1] .- state.height[:,1] .* state.grad_h[:,2])./((state.height[:,1] .+ state.height[:,2]) .* (state.height[:,1] .+ state.height[:,2])))./ state.height[:,1]
    state.vel[:,1] .-= sys.D .* ( (state.height[:,2] .* state.grad_h[:,1] .- state.height[:,1] .* state.grad_h[:,2])./((state.height[:,1] .+ state.height[:,2]) ))./ state.height[:,1]
    # state.vel[:,2] .-= sys.D .* ( (state.height[:,1] .* state.grad_h[:,2] .- state.height[:,2] .* state.grad_h[:,1])./((state.height[:,1] .+ state.height[:,2]) .* (state.height[:,1] .+ state.height[:,2])))./ state.height[:,2]
    state.vel[:,2] .-= sys.D .* ( (state.height[:,1] .* state.grad_h[:,2] .- state.height[:,2] .* state.grad_h[:,1])./((state.height[:,1] .+ state.height[:,2]) ))./ state.height[:,2]
    return nothing
end

"""
	function moments!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)

same as [moments!](@ref) but for multlayer systems
"""
function moments!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)
    # Get views of the populations
    f0, f1, f2 = Swalbe.viewdistsMultiLayer_1D(state.fout) 
    # Compute the height
    for i in 1:sys.layers
        height_view=view(state.height,:,i)
        sum!(height_view, state.fout[:,:,i])
    end
    # and the velocities (as simple as possible)
    state.vel .= (f1 .- f2) ./ state.height
    return nothing
end
