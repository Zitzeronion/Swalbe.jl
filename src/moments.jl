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
function moments!(height::Matrix, velx, vely, fout)
    # Get views of the populations
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(fout) 
    # Compute the height
    sum!(height, fout)
    # and the velocities (as simple as possible)
    velx .= (f1 .- f3 .+ f5 .- f6 .- f7 .+ f8) ./ height
    vely .= (f2 .- f4 .+ f5 .+ f6 .- f7 .- f8) ./ height
    return nothing
end
# with new state struct
function moments!(state::State)
    # Get views of the populations
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    # and the velocities (as simple as possible)
    state.velx .= (f1 .- f3 .+ f5 .- f6 .- f7 .+ f8) ./ state.height
    state.vely .= (f2 .- f4 .+ f5 .+ f6 .- f7 .- f8) ./ state.height
    return nothing
end
#with active thin film model
function moments!(state::StateActive)
    # Get views of the populations
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    # and the velocities (as simple as possible)
    state.velx .= (f1 .- f3 .+ f5 .- f6 .- f7 .+ f8) ./ state.height
    state.vely .= (f2 .- f4 .+ f5 .+ f6 .- f7 .- f8) ./ state.height
    return nothing
end

function moments!(state::CuState)
    # Get views of the populations
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(state.fout) 
    # Compute the height
    # TODO: figuring out this new CUDA problem, seems `sum!` is broken
    sum!(state.height, state.fout)[:,:,1]
    # and the velocities (as simple as possible)
    state.velx .= (f1 .- f3 .+ f5 .- f6 .- f7 .+ f8) ./ state.height
    state.vely .= (f2 .- f4 .+ f5 .+ f6 .- f7 .- f8) ./ state.height
    return nothing
end

function moments!(height, vel, fout)
    # Get views of the populations
    f0, f1, f2 = Swalbe.viewdists_1D(fout) 
    # Compute the height
    sum!(height, fout)
    # and the velocities (as simple as possible)
    vel .= (f1 .- f2) ./ height
    return nothing
end

function moments!(state::State_1D)
    # Get views of the populations
    f0, f1, f2 = Swalbe.viewdists_1D(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    # and the velocities (as simple as possible)
    state.vel .= (f1 .- f2) ./ state.height
    return nothing
end
function moments_new!(state::State_1D)
    # Get views of the populations
    # f0, f1, f2 = Swalbe.viewdists_1D(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    return nothing
end
function moments_overdamped_01!(state::State_1D)
    # Get views of the populations
    # f0, f1, f2 = Swalbe.viewdists_1D(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    state.vel .= -state.h∇p .* state.slip
    return nothing
end
function moments_overdamped_23!(state::State_1D)
    # Get views of the populations
    # f0, f1, f2 = Swalbe.viewdists_1D(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    state.vel .= -state.h∇p ./ state.slip
    return nothing
end
#active thin film version
function moments!(state::Active_1D)
    # Get views of the populations
    f0, f1, f2 = Swalbe.viewdists_1D(state.fout) 
    # Compute the height
    sum!(state.height, state.fout)
    # and the velocities (as simple as possible)
    state.vel .= (f1 .- f2) ./ state.height
    return nothing
end

"""
    rho_moments!(state)
Calculates the moments for ``\\rho`` as 

`` \\rho = \\sum g_i`` 

`` \\rho u = v - D\\frac{\\nabla \\rho}{\\rho}- DA_2 \\nabla \\rho + D \\nabla h \\nabla_h \\mathcal{F}``

That leads to solving 

`` \\partial_t \\rho = \\nabla \\cdot\\left[ \\rho^2 u\\right]

``A_2`` is a second order diffusive term that standarized is choosen ``A_2=0``. Also we iclude a window integration to prevent numerical noise at step functions

`` \\rho(x) \\mapsto \\rho \\omega+ \\frac{1-\\omega}{2}(\\rho(x+\\Delta x) + \\rho(x-\\Delta x)``  

where standarized `` \\omega =0``. 

That is like the standard model for advection difusion as seen in [Tim Krueger Lattice Boltzman methods] but without the diffusive term in the equilibrium function. 

See also: [`Swalbe.rho_equilibrium`](@ref)

"""
function rho_moments!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho, state.gout)
    # Window integration, as standard value `sys.weight_center=0` does that cost computation time?
    # rip, rim = viewneighbors_1D(state.dgrad)
    # circshift!(rip, state.rho, 1)
    # circshift!(rim, state.rho, -1)
    # state.rho .= state.rho  .* sys.weight_center .+ (1-sys.weight_center)/2 .* rip .+ (1-sys.weight_center)/2 .* rim 
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho, state.rho, state.dgrad)
     ∇f!(state.grad_p, state.pressure, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
     state.D .= sys.D .* ( 1 .+ sys.A_2 .* state.rho)
    #  state.D .= sys.D .* ( 1 .+ sys.A_2 .* max.(0, 50 .* state.rho .- 500))
    #  state.D .= sys.D .*(1 .+ sys.A_2 .* ifelse.(state.rho.>10, -200/7 .- 5/7 .* state.rho .+ 5/14 .* state.rho .* state.rho, 0))
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


 function rho_moments_0!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho, state.gout)
    # Window integration, as standard value `sys.weight_center=0` does that cost computation time?
    # rip, rim = viewneighbors_1D(state.dgrad)
    # circshift!(rip, state.rho, 1)
    # circshift!(rim, state.rho, -1)
    # state.rho .= state.rho  .* sys.weight_center .+ (1-sys.weight_center)/2 .* rip .+ (1-sys.weight_center)/2 .* rim 
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho, state.rho, state.dgrad)
     ∇f!(state.grad_p, state.pressure, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
     state.D .= sys.D .* ( 1 .+ sys.A_2 .* state.rho)
    #  state.D .= sys.D .* ( 1 .+ sys.A_2 .* max.(0, 50 .* state.rho .- 500))
    #  state.D .= sys.D .*(1 .+ sys.A_2 .* ifelse.(state.rho.>10, -200/7 .- 5/7 .* state.rho .+ 5/14 .* state.rho .* state.rho, 0))
        state.rho_adv_vel .= (
            state.∇γ ./ sys.μ .* state.V_gamma
            .- state.grad_p ./ sys.μ .* state.V_p
        )
            state.rho_vel .= state.rho_adv_vel  .- state.D .* state.grad_rho ./ state.rho  .- sys.D .* state.grad_h .* state.nabla_F
     return nothing
 end

"""
    rho_A_moments!(state)
Calculates the moments for ``\\rho_A`` as 

`` \\rho_A = \\sum h_i`` 

`` \\rho_A u = v - D_A\\frac{\\nabla \\rho_A}{\\rho_A}+ D_A \\nabla h \\nabla_h \\mathcal{F}``

That leads to solving 

`` \\partial_t \\rho_A = \\nabla \\cdot\\left[ \\rho_A^2 u\\right]


We only have implemented ``\\alpha=0`` for ``\\rho_A`` i.e. we assume the product to be neutrally buoyant. 

That is like the standard model for advection difusion as seen in [Tim Krueger Lattice Boltzman methods] but without the diffusive term in the equilibrium function. 

See also: [`Swalbe.rho_equilibrium`](@ref)

"""
function rho_A_moments!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho_A, state.hout)
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho_A, state.rho_A, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
    #  state.D .= sys.D .* ( 1 .+ sys.A_2 .* state.rho)
    #  state.D .= sys.D .* ( 1 .+ sys.A_2 .* max.(0, 50 .* state.rho .- 500))
    #  state.D .= sys.D .*(1 .+ sys.A_2 .* ifelse.(state.rho.>10, -200/7 .- 5/7 .* state.rho .+ 5/14 .* state.rho .* state.rho, 0))
     state.rho_A_vel .= state.vel .- sys.D_Product .* state.grad_rho_A ./ state.rho_A .+ sys.D_Product .* state.grad_h ./ (state.height)
     return nothing
 end
 

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

 function rho_A_moments_0!(state::Active_1D, sys::SysConstActive_1D)
     # Compute the density
     sum!(state.rho_A, state.hout)
     # calculate derivatives
     ∇f!(state.grad_h, state.height, state.dgrad)
     ∇f!(state.grad_rho_A, state.rho_A, state.dgrad)
     # and the velocities (forcefeeded by the model equations)
        state.rho_A_adv_vel .= (
            state.∇γ ./ sys.μ .* state.V_gamma_A
            .- state.grad_p ./ sys.μ .* state.V_p_A
        )
     state.rho_A_vel .= state.rho_A_adv_vel .- sys.D_Product .* state.grad_rho_A ./ state.rho_A .+ sys.D_Product .* state.grad_h ./ (state.height)
     return nothing
 end
"""
	function rho_moments_linear!(state::Active_1D)


linear verion of the LBM schem for active thin films. Untested and unused. 
"""
function rho_moments_linear!(state::Active_1D)
    # Get views of the populations
    g0, g1, g2 = Swalbe.viewdists_1D(state.gout) 
    # Compute the height
    sum!(state.rho, state.gout)
    # and the velocities (as simple as possible)
    state.rho_vel .= (g1 .- g2) ./ state.rho
    return nothing
end
