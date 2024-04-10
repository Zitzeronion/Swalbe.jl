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
One well studied assumption is that the fluid velocity vanishes at `` h(\\mathbf{x}) = 0 `` which is called **no-slip** condition.
In terms of thin film mobility `M(h)` this means 

`` M(h) = \\frac{h^3}{3\\mu}. ``

This assumption, however, can be relaxed to allow for some slippage with a further parameter δ which determines the slip length.
The here presented model matches the thin film equation by modification of the momentum equation of the shallow water model.
Therefore the force we add to the momentum equation to compensate for the substrate fluid friction is given as

`` F_{fric} = -\\nu \\alpha_{\\delta}(h)\\mathbf{u},  ``

where ν is the kinematic viscosity and α 

`` \\alpha_{\\delta}(h) = \\frac{6h}{2h^2 + 6h\\delta +3\\delta^2}, ``

and δ being the slip length, or the distance inside the substrate where the fluids velocity vanishes.

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
# with state struct
function slippage!(state::State, sys::SysConst)
    state.slipx .= (6*sys.μ .* state.height .* state.velx) ./ (2 .* state.height .* state.height .+ 6sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    state.slipy .= (6*sys.μ .* state.height .* state.vely) ./ (2 .* state.height .* state.height .+ 6sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    return nothing
end
#for active thin film
function slippage!(state::StateActive, sys::SysConstActive)
    state.slipx .= (6*sys.μ .* state.height .* state.velx) ./ (2 .* state.height .* state.height .+ 6sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    state.slipy .= (6*sys.μ .* state.height .* state.vely) ./ (2 .* state.height .* state.height .+ 6sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    return nothing
end

function slippage!(slip, height, vel, δ, μ)
    slip .= (6μ .* height .* vel) ./ (2 .* height.^2 .+ 6δ .* height .+ 3δ^2 )
    return nothing
end
# with state struct
function slippage!(state::State_1D, sys::SysConst_1D)
    state.slip .= (6*sys.μ .* state.height .* state.vel) ./ (2 .* state.height .* state.height .+ 6*sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    return nothing
end
#for active thin film
function slippage!(state::Active_1D, sys::SysConstActive_1D)
    state.slip .= (6*sys.μ .* state.height .* state.vel) ./ (2 .* state.height .* state.height .+ 6*sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    return nothing
end
# Seems bullshit
function slippage_new!(state::State_1D, sys::SysConst_1D)
    state.slip .= (6*sys.μ)./(2 .* state.height.^3 .+ 6*sys.δ .* state.height .* state.height .+ 3*sys.δ*sys.δ .* state.height)
    return nothing
end

"""
    h∇p!(state)

Computation of the pressure gradient multiplied with the height.

# Mathematics

The pressure gradient (``\\nabla p``) is **the** driving force of the standard thin film equation `` \\partial_t h = (M(h)\\nabla p)``.
Our approach however does not solve the thin film equation directly,
Instead we have to add the pressure gradient as a force which is given as

`` F_{film} = -\\frac{1}{\\rho_0} h \\nabla p_{film}, ``

where the term ``p_{film}`` describes the film pressure 

`` p_{film} = -\\gamma [\\Delta h - \\Pi(h)] .``

As such it is the combination of a laplacian term that minimizes the surface area as well as a interfacial potential between substrate and fluid.

# Examples
```jldoctest
julia> using Swalbe, Test

julia> state = Swalbe.Sys(Swalbe.SysConst(Lx=10, Ly=10), "CPU");

julia> state.pressure .= reshape(collect(1:100),10,10);

julia> Swalbe.h∇p!(state)

julia> @test all(state.h∇px[1,:] .== -4.0) # at boundary
Test Passed

julia> @test all(state.h∇px[2,:] .== 1.0) # inside
Test Passed

```

# References

- [Zitz, Scagliarini and Harting](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)
- [Craster, Matar](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131)
- [Oron, Davis and Bankoff](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.69.931)

See also: [`Swalbe.filmpressure!`](@ref)
"""
function h∇p!(state::State)
    fip, fjp, fim, fjm, fipjp, fimjp, fimjm, fipjm = viewneighbors(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(fip, state.pressure, (1,0))
    circshift!(fjp, state.pressure, (0,1))
    circshift!(fim, state.pressure, (-1,0))
    circshift!(fjm, state.pressure, (0,-1))
    # Diagonal elements  
    circshift!(fipjp, state.pressure, (1,1))
    circshift!(fimjp, state.pressure, (-1,1))
    circshift!(fimjm, state.pressure, (-1,-1))
    circshift!(fipjm, state.pressure, (1,-1))
    # In the end it is just a weighted sum...
    state.h∇px .= state.height .* (-1/3 .* (fip .- fim) .- 1/12 .* (fipjp .- fimjp .- fimjm .+ fipjm))
    state.h∇py .= state.height .* (-1/3 .* (fjp .- fjm) .- 1/12 .* (fipjp .+ fimjp .- fimjm .- fipjm))

    return nothing
end
#Active thin film version
function h∇p!(state::StateActive)
    fip, fjp, fim, fjm, fipjp, fimjp, fimjm, fipjm = viewneighbors(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(fip, state.pressure, (1,0))
    circshift!(fjp, state.pressure, (0,1))
    circshift!(fim, state.pressure, (-1,0))
    circshift!(fjm, state.pressure, (0,-1))
    # Diagonal elements  
    circshift!(fipjp, state.pressure, (1,1))
    circshift!(fimjp, state.pressure, (-1,1))
    circshift!(fimjm, state.pressure, (-1,-1))
    circshift!(fipjm, state.pressure, (1,-1))
    # In the end it is just a weighted sum...
    state.h∇px .= state.height .* (-1/3 .* (fip .- fim) .- 1/12 .* (fipjp .- fimjp .- fimjm .+ fipjm))
    state.h∇py .= state.height .* (-1/3 .* (fjp .- fjm) .- 1/12 .* (fipjp .+ fimjp .- fimjm .- fipjm))

    return nothing
end

# One dimensional implementation
function h∇p!(state::State_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    
    # In the end it is just a weighted sum...
    state.h∇p .= state.height .* -0.5 .* (fip .- fim)

    return nothing
end
function h∇p_new!(state::State_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    
    # In the end it is just a weighted sum...
    state.h∇p .= -0.5 .* (fip .- fim)

    return nothing
end
#active thin film version
function h∇p!(state::Active_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    
    # In the end it is just a weighted sum...
    state.h∇p .= state.height .* -0.5 .* (fip .- fim)

    return nothing
end

"""
	function rho_grad_p!(state::Active_1D)

Calculates ``\\rho\\nabla p`` where p is the pressure for ``\\rho``

See also: [`Swalbe.rho_grad_p!`](@ref)
"""
function rho_grad_p!(state::Active_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    
    # In the end it is just a weighted sum...
    state.rho∇p .= state.rho .* -0.5 .* (fip .- fim)
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

The classical thin film equation is an equation without thermal noise which is defined as

``\\partial_t h = \\nabla\\cdot(M(h)\\nabla p).``

Classically thermal excitations are neglected for thin film flows.
It is of course possible to introduce for example an envolving surface tension or viscosity but both do not account for e.g. the spectrum of thermal capillary waves.
The reason for this short coming is the complex from the thin film equation takes if derived from the Landau-Lifshitz Navier-Stokes equation.
One addition term due to the stochastic stress tensor makes it somewhat impossible to solve the eqaution self-consistently.
However Grün et al. showed that it is not nessary to use the full stochastic stress tensor, but simple a multiplicative noise term

``\\partial_t h = \\nabla\\cdot[M(h)(\\nabla p + \\sigma \\mathcal{N}], ``

where ``\\sigma`` and ``\\mathcal{N}`` are a dimensionless temperture and a gaussian white noise.
Similar to the pressure gradient the addition of this term is introduced as force in our model

`` F_{fluc} = \\frac{1}{\\rho_0}\\sqrt{2k_BT\\mu\\alpha_{\\delta}(h)}\\mathcal{N},``

with ``\\alpha_{\\delta}(h)`` being the force generated due to substrate slip. 

# Examples
```jldoctest
julia> using Swalbe, Statistics, Test

julia> x = ones(50,50); y = ones(50,50); h = ones(50,50);

julia> Swalbe.thermal!(x, y, h, 0.1, 1/6, 1)

julia> @test mean(x) ≈ 0.0 atol=1e-2
Test Passed

julia> @test mean(y) ≈ 0.0 atol=1e-2
Test Passed

julia> @test var(x) ≈ 2*0.1/11 atol=(2*0.1/11)/10 # var = 2kbt*6*μ/slip
Test Passed

```

# References

- [Grün, Mecke and Rauscher](https://link.springer.com/article/10.1007/s10955-006-9028-8)
- [Mecke, Rauscher](https://iopscience.iop.org/article/10.1088/0953-8984/17/45/042/meta)
- [Davidovitch, Moro and Stone](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.244505)

See also: [`Swalbe.slippage!`](@ref)
""" 
function thermal!(fluc_x, fluc_y, height, kᵦT, μ, δ)
    randn!(fluc_x)
    randn!(fluc_y)
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

function thermal!(fluc, height, kᵦT, μ, δ)
    randn!(fluc)
    fluc .*= sqrt.(2 .* kᵦT .* μ .* 6 .* height ./
                  (2 .* height.^2 .+
                   6 .* height .* δ .+
                   3 .* δ^2))
    
    return nothing
end

"""
   inclination!(α, state)

Force that mimics the effect of an inclined plate.

Simple model to add a body force on the fluid that should mimic the effect of an inclined plate.
The model includes a smoothing `tanh` function to absorb shocks which might occure.
    
# Arguments

-`α :: Vector`: Force vector that both hosts direction of the force as well as strength
-`state::State`: Lattice Boltzmann state of the fluid, here we need the `state.Fx`, `state.Fy` fields
-`t::Int`: Time step, used for the smoothing `tanh` factor
-`tstart::Int`: Time delay at which the `tanh` becomes positve
-`tsmooth::Int`: Time interval over which the `tanh` is smeared

# Mathematics

This body force is simply force strength time the mass of the fluid or even simpler the height (assuming ρ=1).
Thus the force becomes

`` \\mathbf{F} = \\mathbf{\\alpha} h \\tanh\\bigg(\\frac{t-t_0}{t_s}\\bigg),``

with ``t_0`` being the time lag at which the `tanh` changes sign and ``t_s`` is width of interval between -1 and 1.

See also: [Swalbe.run_dropletforced](@ref)
"""
function inclination!(α::Vector, state::State; t=1000, tstart=0, tsmooth=1)
    state.Fx .+= state.height .* α[1] .* (0.5 .+ 0.5 .* tanh((t - tstart)/tsmooth))
    state.Fy .+= state.height .* α[2] .* (0.5 .+ 0.5 .* tanh((t - tstart)/tsmooth))

    return nothing
end

function inclination!(α::Vector, state::State_1D; t=1000, tstart=0, tsmooth=1)
    state.F .+= state.height .* α[1] .* (0.5 .+ 0.5 .* tanh((t - tstart)/tsmooth))

    return nothing
end



"""
    surface_tension_gradient!(state)

Computes the gradient of a spatially resolved surface tension field via central differneces and then concludes the surface stress as 

`` \\frac{h^2 + 2 \\delta h }{2M(h)}``

with 

``M(h)= \\frac{2h^2 + 6 \\delta h + 3\\delta^2 }{ 6}``

and implicitly ``\\rho_h=1``. 
"""
function ∇γ!(state::Active_1D, sys::SysConstActive_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.γ, 1)
    circshift!(fim, state.γ, -1)
    
    # In the end it is just a weighted sum...
    # state.∇γ .= -3/2 .* ((fip .- fim) ./ 2.0)
    state.∇γ .=((fip .- fim) .* -0.5)
    state.stress .= state.∇γ .* (state.height .* state.height .* 0.5 .+ sys.δ .* state.height) .* 6 ./ (2 .* state.height .* state.height .+ 6*sys.δ .* state.height .+ 3*sys.δ*sys.δ )
    return nothing
end





"""
    rho_forcing!(state,sys)
Calculates the neccessary forces for `` \\rho `` according to 

`` \\F = - \\frac{D \\nabla \\rho + M\\Gamma \\rho \\nabla \\rho - D \\rho \\nabla \\ln h}{\\rho + \\rho^*}``

where ``\\rho^*`` ensure we don't devide by ``0``.

All unused. Maybe delete them.
"""

function rho_forcing!(state::CuStateActive_1D,sys)
     # calculate derivatives as finite differences
     ∇f_Cu_1D!(state.grad_rho, state.rho, state.dgrad)
     ∇f_Cu_1D!(state.grad_h, state.height, state.dgrad)
     # calculate force
     state.rho_F .= ( (-sys.D .* state.grad_rho .- sys.M .* sys.Γ .* state.rho .* state.grad_rho .+ sys.D .* state.rho .* (state.grad_h ./ (state.height .+ 2*sys.hcrit))) ./ (state.rho .+ 2*sys.hcrit))
    return nothing
end

function rho_forcing!(state::StateActive_1D,sys)
    # calculate derivatives as finite differences
    ∇f!(state.grad_rho, state.rho, state.dgrad)
    ∇f!(state.grad_h, state.height, state.dgrad)
    # calculate force
    state.rho_F .= - sys.D .* state.grad_rho .+state.vel .* state.rho .+ sys.D .* state.rho .* state.grad_h ./ state.height .- state.rho.*state.rho_vel
   return nothing
end

function rho_diffusion!(state::StateActive_1D,sys)
    # calculate derivatives as finite differences
    ∇f!(state.grad_rho, state.rho, state.dgrad)
    ∇f!(state.grad_h, state.height, state.dgrad)
    # calculate force
    state.diffusion_rho .=  sys.D .* state.grad_rho .- sys.D .* state.rho .* state.grad_h ./ state.height
   return nothing
end


"""
    rho_forcing!(state,sys)
Calculates the neccessary forces for `` \\rho `` according to 

`` F = \\rho \\nabla p(\\rho) ``

where ``p(\\rho)=(0.05/(0.05+\\rho))^n`` .
"""
function rho_forcing_quadratic!(state::StateActive_1D,sys)
    # add a force coming from a disjoining pressure p(rho)=0.05^15/(rho+0.05)^15 as F= rho nabla p(rho)
    # carfull there is another call of this line line rho_moments_quadratic!
    ∇f!(state.grad_rho, state.rho, state.dgrad)
    state.rho_F .= sys.rho_n .* Swalbe.power_broad.(sys.hcrit./(state.rho .+ sys.hcrit), sys.rho_n) ./ (state.rho .+ sys.hcrit) .* state.grad_rho
   return nothing
end

"""
    function rho_grad_disj_p!(state::Active_1D)

A central differences scheme to calculate ``\\rho \\nabla p_{\\rho}``.
"""
function rho_grad_disj_p!(state::Active_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.rho_pressure, 1)
    circshift!(fim, state.rho_pressure, -1)
    
    # In the end it is just a weighted sum...
    state.rho∇p .= state.rho .* 0.5 .* (fip .- fim)

    return nothing
end


"""
    function rho_A_grad_disj_p!(state::Active_1D)

A central differences scheme to calculate ``\\rho_A \\nabla p_{\\rho}``.
"""
function rho_A_grad_disj_p!(state::Active_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.rho_A_pressure, 1)
    circshift!(fim, state.rho_A_pressure, -1)
    
    # In the end it is just a weighted sum...
    state.rho_A∇p .= state.rho_A .* 0.5 .* (fip .- fim)

    return nothing
end



"""
    function rho_B_grad_disj_p!(state::Active_1D)

A central differences scheme to calculate ``\\rho_B \\nabla p_{\\rho}``.
"""
function rho_B_grad_disj_p!(state::Active_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.rho_B_pressure, 1)
    circshift!(fim, state.rho_B_pressure, -1)
    
    # In the end it is just a weighted sum...
    state.rho_B∇p .= state.rho_B .* 0.5 .* (fip .- fim)

    return nothing
end



function grad_rho_u!(state::Active_1D, sys::SysConstActive_1D)
    state.z.=max.(state.height .- sys.B, 0)
    state.rho_adv_vel .= state.rho .* ( state.h∇p ./ (state.height.*sys.μ)  .* ( state.z.^2 ./ 2 .- state.height .* state.z .- sys.δ .* state.height .- sys.δ ./2) .+ state.∇γ ./ sys.μ .* (state.z .+ sys.δ))
    # state.rho_adv_vel .= state.rho .* state.vel
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.rho_adv_vel, 1)
    circshift!(fim, state.rho_adv_vel, -1)
    
    # In the end it is just a weighted sum...
    state.∇rho_u .= 0.5 .* (fip .- fim)
    # state.∇rho_u .= 0

    return nothing
end




function rho_grad_v!(state::Active_1D, sys::SysConstActive_1D)
    state.rho_adv_vel .= state.vel .- sys.D .* state.grad_rho ./ state.rho .+ sys.D .* state.grad_h ./ (state.height)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.rho_adv_vel, 1)
    circshift!(fim, state.rho_adv_vel, -1)
    state.diffusion_rho.= state.rho.*(fip .- fim)
    return nothing
end



function rho_forcing_new_old_approach!(state::Active_1D, sys::SysConstActive_1D)
	∇f!(state.grad_rho, state.rho, state.dgrad)
	∇f!(state.grad_h, state.height, state.dgrad)
	∇f!(state.grad_p, state.pressure, state.dgrad)
	state.rho_F .= (
		.- state.rho .* state.grad_p 
		.+ ( .+ sys.Γ .* state.rho .* state.grad_rho .- sys.D .* state.grad_rho .* sys.μ .+ sys.D .* state.rho .* state.grad_h ./ state.height .* sys.μ) .* (
				(sys.δ .+ state.height ./ 2) ./ (1/6 .* (3 * sys.δ*sys.δ .+ 6 .* sys.δ .* state.height .+ 2 .* state.height .* state.height))
		) 
		.- state.rho .* state.rho_vel .* sys.μ ./  (1/6 .* (3 * sys.δ*sys.δ .+ 6 .* sys.δ .* state.height .+ 2 .* state.height .* state.height))
	)
    return nothing
end



function rho_forcing_new_approach!(state::Active_1D, sys::SysConstActive_1D)
        ∇f!(state.grad_rho, state.rho, state.dgrad)
        ∇f!(state.grad_h, state.height, state.dgrad)
        ∇f!(state.grad_p, state.pressure, state.dgrad)
        state.rho_F .= (
                - state.rho .* state.grad_p
                .+ ( - sys.Γ .* state.rho .* state.grad_rho  .- sys.D .* state.grad_rho .* sys.μ  .- sys.D .* state.rho .* state.grad_h .* state.nabla_F .* sys.μ ) .* ( state.V_gamma ./ state.V_p)
                .- state.rho .* state.rho_vel .* sys.μ ./  ( state.V_p )
	)
    return nothing
end



function rho_forcing_new_approach_bounce_back!(state::Active_1D, sys::SysConstActive_1D)
        ∇f!(state.grad_rho, state.rho, state.dgrad)
        ∇f!(state.grad_h, state.height, state.dgrad)
        ∇f!(state.grad_p, state.pressure, state.dgrad)
        state.rho_F .= (
                - state.rho .* state.grad_p
                .+ ( - sys.Γ .* state.rho .* state.grad_rho  .- sys.D .* state.grad_rho .* sys.μ  .- sys.D .* state.rho .* state.grad_h .* state.nabla_F .* sys.μ ) .* ( state.V_gamma ./ state.V_p)
                .- state.rho .* state.rho_vel .* sys.μ ./  ( state.V_p )
		.+ state.rho∇p
        )
        state.rho_F .*= (1 .- state.rightrocks) .* (1 .- state.leftrocks)
    return nothing
end

"""
    function V_gamma!(state::Active_1D, sys::SysConstActive_1D)

Calculates the term ``V_\\gamma`` that enters as a part of the diffusive velocity of the solvent ``\\rho`` as ``V_\\gamma \\nabla \\gamma``. We have 

``V_\\gamma = \\begin{cases}
\\delta & \\alpha = \\infty\\\\
\\frac{1}{1-e^{-\\alpha h}}\\left( -\\delta e^{-\\alpha h}+ \\frac{1-e^{-\\alpha h }( \\alpha h +1 )}{\\alpha}+\\delta\\right) & \\alpha \\geq 0 \\\\
\\delta + \\frac{h}{2} & \\alpha =0\\\\
\\frac{1}{e^{\\alpha h}-1}\\left( \\delta (e^{\\alpha h}-1)+ \\frac{e^{\\alpha h}-1}{\\alpha }-h \\right) & \\alpha \\leq 0\\\\
\\delta + h & \\alpha =-\\infty
\\end{cases}``
"""
function V_gamma!(state::Active_1D, sys::SysConstActive_1D)
    if sys.alpha==-Inf
        state.V_gamma .= sys.δ .+ state.height
    elseif sys.alpha<0
        state.V_gamma .= (sys.δ .* (exp.(sys.alpha .* state.height) .- 1) .+ (exp.(sys.alpha .* state.height) .- 1) ./ sys.alpha .- state.height ) ./ (exp.(sys.alpha .* state.height) .- 1)
    elseif sys.alpha==0
        state.V_gamma .= (2 .* sys.δ .+ state.height ) ./ 2
    elseif sys.alpha>0 && sys.alpha < Inf
        state.V_gamma .=  ( -sys.δ .* exp.(-sys.alpha .* state.height) .+ (1 .- exp.(-sys.alpha .* state.height) .* (sys.alpha .* state.height .+ 1)) ./ sys.alpha .+ sys.δ) ./ (1 .- exp.(-sys.alpha .* state.height))
    else
        state.V_gamma .=  sys.δ
    end
    return nothing
end

function V_gamma_A!(state::Active_1D, sys::SysConstActive_1D)
        state.V_gamma_A .= (2 .* sys.δ .+ state.height ) ./ 2
    return nothing
end




"""
    function V_p!(state::Active_1D, sys::SysConstActive_1D)

Calculates ``V_p`` that determines the pressure driven part of the solvent advection that has a velocity ``\\frac{\\nalbla p }{\\mu }V_p`` as

``V_p= \\begin{cases}
\\frac{\\delta^2}{2}+ \\delta h & \\alpha =\\infty\\\\
- \\frac{2 - 2 \\alpha h - \\alpha^2 \\delta (\\delta  + 2 h) + e^{-\\alpha h} ( \\alpha^2 (\\delta  + h)^2-2)}{2 \\alpha^2(1-e^{-\\alpha h})} & \\alpha \\geq 0 \\\\
\\frac{1}{6} (3 \\delta^2 + 6 \\delta h + 2 h^2) & \\alpha =0 \\\\
-  \\frac{-2 + \\alpha^2 (\\delta + h)^2 - e^{\\alpha h}\\left(-2 + 2 \\alpha h + \\alpha^2 \\delta (\\delta + 2 h)\\right)}{2 a^2(e^{\\alpha h}-1)} & \\alpha \\leq 0\\
\\frac{1}{2}(\\delta^2 + 2\\delta + h^2) & \\alpha =-\\infty
\\end{cases}``
"""
function V_p!(state::Active_1D, sys::SysConstActive_1D)
    if sys.alpha==-Inf
        state.V_p .= (sys.δ*sys.δ .+ 2 * sys.δ .* state.height .+ state.height .* state.height) .* 0.5
    elseif sys.alpha<0
        state.V_p .= .-(-2 .+ sys.alpha*sys.alpha .* (sys.δ .+ state.height) .* (sys.δ .+ state.height).- exp.(sys.alpha .* state.height) .* (-2 .+ 2*sys.alpha .* state.height .+ sys.alpha*sys.alpha * sys.δ .*(sys.δ .+ 2 .* state.height))) ./ (2*sys.alpha*sys.alpha .* (exp.(sys.alpha .* state.height) .- 1))
    elseif sys.alpha==0
        state.V_p .= (2 .* state.height .* state.height .+ 6*sys.δ .* state.height .+ 3*sys.δ*sys.δ ) ./ (6)
    elseif sys.alpha>0 && sys.alpha < Inf
        state.V_p .= .-(2 .-2 * sys.alpha .* state.height .- sys.alpha*sys.alpha * sys.δ .* (sys.δ .+ 2 .* state.height ) .+ exp.(-sys.alpha .* state.height) .* (sys.alpha*sys.alpha .*(sys.δ .+ state.height) .* (sys.δ .+ state.height).-2)) ./ (2*sys.alpha*sys.alpha .* (1 .- exp.(-sys.alpha .* state.height)))
    else
        state.V_p .= sys.δ*sys.δ/2  .+ sys.δ .* state.height
    end
    return nothing
end

function V_p_A!(state::Active_1D, sys::SysConstActive_1D)
        state.V_p_A .= (2 .* state.height .* state.height .+ 6*sys.δ .* state.height .+ 3*sys.δ*sys.δ ) ./ (6)
    return nothing
end



"""
    function nabla_F!(state::Active_1D, sys::SysConstActive_1D)

Calculates the Fick Jackobs term that corrects the diffusion according to a exponential distribution of solute. 

``\\beta \\frac{1}{\\nabla h} \\nabla \\mathcal{F}= \\begin{cases}
0 & \\alpha = \\pm \\infty\\\\
- \\frac{1}{h} & \\alpha =0\\\\
\\frac{\\alpha}{e^{\\alpha h}-1} & \\alpha \\geq 0\\\\
\\frac{\\alpha e^{\\alpha h}}{e^{\\alpha h}-1} & \\alpha \\leq 0  
\\end{cases}``
"""
function nabla_F!(state::Active_1D, sys::SysConstActive_1D)
    if sys.alpha==-Inf
        state.nabla_F .= 0
    elseif sys.alpha<0
        state.nabla_F .=  .- sys.alpha .* exp.(sys.alpha .* state.height) ./ (exp.(sys.alpha .* state.height) .- 1)
    elseif sys.alpha==0
        state.nabla_F .= .-1 ./ state.height
    elseif sys.alpha>0 && sys.alpha < Inf
        state.nabla_F .=  .- sys.alpha ./ (exp.(sys.alpha .* state.height) .- 1)
    else
        state.nabla_F .= 0
    end
    return nothing
end



function nabla_F_precursor!(state::Active_1D, sys::SysConstActive_1D)
    if sys.alpha==-Inf
        state.nabla_F .= 0
    elseif sys.alpha<0
        state.nabla_F .=  .- sys.alpha .* exp.(sys.alpha .* state.height) ./ (exp.(sys.alpha .* state.height) .- 1)
    elseif sys.alpha==0
        state.nabla_F .= .-1 ./ state.height
    elseif sys.alpha>0 && sys.alpha < Inf
        state.nabla_F .=  .- sys.alpha ./ (exp.(sys.alpha .* state.height) .- 1)
    else
        state.nabla_F .= 0
    end
    @. state.precursor = ifelse(state.height <= sys.hmin - sys.hcrit + 0.01, 0,1)
    state.nabla_F .-= (1 .- state.precursor) .* sys.prescursor_nabla_F
    return nothing
end

"""
    function nabla_F_0!(state::Active_1D, sys::SysConstActive_1D)

Calculates the Fick Jackobs term that corrects the diffusion according to a exponential distribution of solute. 

``\\beta \\frac{1}{\\nabla h} \\nabla \\mathcal{F}= - \\frac{1}{h} ``
"""
function nabla_F_0!(state::Active_1D)
    state.nabla_F_0 .= .-1 ./ state.height
    return nothing
end



"""
    function V_p_0!(state::Active_1D, sys::SysConstActive_1D)

Calculates ``V_p`` that determines the pressure driven part of the solvent advection that has a velocity ``\\frac{\\nalbla p }{\\mu }V_p`` as

``V_p= \\frac{1}{6} (3 \\delta^2 + 6 \\delta h + 2 h^2)``
"""
function V_p_0!(state::Active_1D, sys::SysConstActive_1D)
    state.V_p_0 .= (2 .* state.height .* state.height .+ 6*sys.δ .* state.height .+ 3*sys.δ*sys.δ ) ./ (6)
    return nothing
end



"""
    function V_gamma_0!(state::Active_1D, sys::SysConstActive_1D)

Calculates the term ``V_\\gamma`` that enters as a part of the diffusive velocity of the solvent ``\\rho_A`` as ``V_\\gamma \\nabla \\gamma``. We have 

``V_\\gamma =\\delta + \\frac{h}{2} & \\alpha =0\\\\``
"""
function V_gamma_0!(state::Active_1D, sys::SysConstActive_1D)
    state.V_gamma .= (2 .* sys.δ .+ state.height ) ./ 2
    return nothing
end
