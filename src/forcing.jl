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
    @. slipx .= (6μ * height * velx) / (2 * height^2 + 6δ * height + 3δ^2 )
    @. slipy .= (6μ * height * vely) / (2 * height^2 + 6δ * height + 3δ^2 )
    return nothing
end

slippage!(state::LBM_state_2D, sys::SysConst) = slippage!(state.slipx, state.slipy, state.height, state.velx, state.vely, sys.param.δ, sys.param.μ)

slippage!(state::Expanded_2D, sys::SysConst) = slippage!(state.basestate.slipx, state.basestate.slipy, state.basestate.height, state.basestate.velx, state.basestate.vely, sys.param.δ, sys.param.μ)

function slippage!(slip, height, vel, δ, μ)
    @. slip .= (6μ * height * vel) / (2 * height^2 + 6δ * height + 3δ^2 )
    return nothing
end

slippage!(state::LBM_state_1D, sys::Consts_1D) = slippage!(state.slip, state.height, state.vel, sys.param.δ, sys.param.μ)
    
slippage!(state::Expanded_1D, sys::Consts_1D) = slippage!(state.basestate.slip, state.basestate.height, state.basestate.vel, sys.param.δ, sys.param.μ)

# Dirty hack for reducing slip length
function slippage2!(state::LBM_state_2D, sys::SysConst)
    @. state.slipx .= (6sys.param.μ * (state.height+sys.param.hcrit) * state.velx) / (2 * (state.height+sys.param.hcrit)^2 + 6sys.param.δ * (state.height+sys.param.hcrit) + 3sys.param.δ^2 )
    @. state.slipy .= (6sys.param.μ * (state.height+sys.param.hcrit) * state.vely) / (2 * (state.height+sys.param.hcrit)^2 + 6sys.param.δ * (state.height+sys.param.hcrit) + 3sys.param.δ^2 )
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

julia> state = Swalbe.Sys(Swalbe.SysConst(Lx=10, Ly=10, param=Swalbe.Taumucs()), "CPU");

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
function h∇p!(state::LBM_state_2D)
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
    @. state.h∇px .= state.height * (-1/3 * (fip - fim) - 1/12 * (fipjp - fimjp - fimjm + fipjm))
    @. state.h∇py .= state.height * (-1/3 * (fjp - fjm) - 1/12 * (fipjp + fimjp - fimjm - fipjm))

    return nothing
end

function h∇p!(state::LBM_state_1D)
    fip, fim = viewneighbors_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    # In the end it is just a weighted sum...
    state.h∇p .= state.height .* -0.5 .* (fip .- fim)

    return nothing
end

function h∇p!(state::Expanded_2D)
    fip, fjp, fim, fjm, fipjp, fimjp, fimjm, fipjm = viewneighbors(state.basestate.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(fip, state.basestate.pressure, (1,0))
    circshift!(fjp, state.basestate.pressure, (0,1))
    circshift!(fim, state.basestate.pressure, (-1,0))
    circshift!(fjm, state.basestate.pressure, (0,-1))
    # Diagonal elements  
    circshift!(fipjp, state.basestate.pressure, (1,1))
    circshift!(fimjp, state.basestate.pressure, (-1,1))
    circshift!(fimjm, state.basestate.pressure, (-1,-1))
    circshift!(fipjm, state.basestate.pressure, (1,-1))
    # In the end it is just a weighted sum...
    @. state.basestate.h∇px .= state.basestate.height * (-1/3 * (fip - fim) - 1/12 * (fipjp - fimjp - fimjm + fipjm))
    @. state.basestate.h∇py .= state.basestate.height * (-1/3 * (fjp - fjm) - 1/12 * (fipjp + fimjp - fimjm - fipjm))

    return nothing
end
# One dimensional implementation
function h∇p!(state::Expanded_1D)
    fip, fim = viewneighbors_1D(state.basestate.dgrad)
    # One dim case, central differences
    circshift!(fip, state.basestate.pressure, 1)
    circshift!(fim, state.basestate.pressure, -1)
    # In the end it is just a weighted sum...
    state.basestate.h∇p .= state.basestate.height .* -0.5 .* (fip .- fim)

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
```jldoctest; output = false
julia> using Swalbe, Statistics, Test, Random

julia> Random.seed!(1234); # Set a seed

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

thermal!(state::State_thermal, sys::SysConst) = thermal!(state.kbtx, state.kbty, state.basestate.height, sys.param.kbt, sys.param.μ, sys.param.δ)

function thermal!(fluc, height, kᵦT, μ, δ)
    randn!(fluc)
    fluc .*= sqrt.(2 .* kᵦT .* μ .* 6 .* height ./
                  (2 .* height.^2 .+
                   6 .* height .* δ .+
                   3 .* δ^2))
    
    return nothing
end

thermal!(state::State_thermal_1D, sys::Consts_1D) = thermal!(state.kbt, state.basestate.height, sys.param.kbt, sys.param.μ, sys.param.δ)

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
function inclination!(α::Vector, state::LBM_state_2D; t=1000, tstart=0, tsmooth=1)
    @. state.Fx .+= state.height * α[1] * (0.5 + 0.5 * tanh((t - tstart)/tsmooth))
    @. state.Fy .+= state.height * α[2] * (0.5 + 0.5 * tanh((t - tstart)/tsmooth))

    return nothing
end

function inclination!(α::Vector, state::Expanded_2D; t=1000, tstart=0, tsmooth=1)
    @. state.basestate.Fx .+= state.basestate.height * α[1] * (0.5 + 0.5 * tanh((t - tstart)/tsmooth))
    @. state.basestate.Fy .+= state.basestate.height * α[2] * (0.5 + 0.5 * tanh((t - tstart)/tsmooth))

    return nothing
end

function inclination!(α::Float64, state::State_1D; t=1000, tstart=0, tsmooth=1)
    state.F .+= state.height .* α .* (0.5 .+ 0.5 .* tanh((t - tstart)/tsmooth))

    return nothing
end

function inclination!(α::Float64, state::Expanded_1D; t=1000, tstart=0, tsmooth=1)
    state.basestate.F .+= state.basestate.height .* α .* (0.5 .+ 0.5 .* tanh((t - tstart)/tsmooth))

    return nothing
end
"""
    update_rho()

Time evolution of the `active` field rho.

TODO: @Tilman!
"""
function update_rho!(rho, rho_int, height, dgrad, differentials; D=1.0, M=0.0)
    lap_rho, grad_rho, lap_h, grad_h = view_four(differentials)
    ∇²f!(lap_rho, rho, dgrad)
    ∇f!(grad_rho, rho, dgrad)
    ∇²f!(lap_h, height, dgrad)
    ∇f!(grad_h, height, dgrad)
    rho_int .= (D .* lap_rho .- 
               M .* (grad_rho.^2 .+ rho .* lap_rho) .- 
               D .* (grad_rho .* grad_h ./ height + rho .* (lap_h ./ height .- (grad_h ./ height).^2 )))

    rho .+= rho_int

    return nothing
end

"""
    surface_tension_gradient!(state)

Computes the gradient of a spatially resolved surface tension field.
"""
function ∇γ!(state::T) where {T <: Expanded_1D}
    fip, fim = viewneighbors_1D(state.basestate.dgrad)
    # One dim case, central differences
    circshift!(fip, state.γ, 1)
    circshift!(fim, state.γ, -1)
    
    # In the end it is just a weighted sum...
    state.∇γ .= -3/2 .* ((fip .- fim) ./ 2.0)
    return nothing
end

"""
    view_four()

Splits a chuck of memory in four equivalent chucks
"""
function view_four(f)
    f1 = view(f, :, 1)
    f2 = view(f, :, 2)
    f3 = view(f, :, 3)
    f4 = view(f, :, 4)
    
    return f1, f2, f3, f4
end
