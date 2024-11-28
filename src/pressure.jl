"""
    filmpressure!(pressure, height, γ, θ, n, m, hmin, hcrit)

Calculation of the capillary pressure which is given by `` p = - γ∇²h+ Π(h) ``.

# Arguments

- `pressure :: Array{Number,2}`: Array that store the result of the compuation
- `height :: Array{Number,2}`: Height field ``h(\\mathbf{x},t)``
- `γ <: Number`: Forcing strenght due to surface tension
- `θ <: Number`: Equilibrium contact angle
- `n :: Int`: Larger power law exponent for `` Π(h) ``
- `m :: Int`: Smaller power law exponent for `` Π(h) ``
- `hmin <: Number`: Parameter of `` Π(h) ``, in fact `` Π(hmin) = 0 ``
- `hcrit <: Number`: Numerical stabilizer for case `` h(\\mathbf{x},t) \\ll hmin ``

# Mathematics

The capillary pressure ``p_{\\text{cap}}`` is the centeral angle to match our model with the thin film equation.
It consists of two parts, first being the laplace pressure `` \\nabla^2 h `` and second being the derivative of the disjoining pontential `` \\Pi(h) ``

`` p_{\\text{cap}} = -\\gamma \\nabla^2 h + \\Pi(h). ``

For the laplacian term we use the same nine point discretization as in `Swlabe.∇²f!`.
`` \\Pi(h) `` on the other hand is given by 

`` \\Pi(h) = \\kappa(\\theta)f(h), ``

where `` \\kappa(\\theta) `` is simply a measure for the **Hamaker constant** and given as

`` \\kappa(\\theta) = \\gamma(1- \\cos(\\theta))\\frac{(n-1)(m-1)}{(n-m)h_{\\text{min}}}.``

For `` f(h) `` one can use various forms, a very common however is the power law given by 

`` f(h) = \\bigg[\\bigg(\\frac{h_{\\text{min}}}{h}\\bigg)^n - \\bigg(\\frac{h_{\\text{min}}}{h}\\bigg)^m\\bigg]. ``

# Examples

```jldoctest
julia> using Swalbe, Test

julia> h = reshape(collect(1.0:25.0),5,5) # A dummy height field
5×5 Matrix{Float64}:
 1.0   6.0  11.0  16.0  21.0
 2.0   7.0  12.0  17.0  22.0
 3.0   8.0  13.0  18.0  23.0
 4.0   9.0  14.0  19.0  24.0
 5.0  10.0  15.0  20.0  25.0

julia> pressure = zeros(5,5); θ = 0.0; # Fully wetting substrate

julia> Swalbe.filmpressure!(pressure, h, zeros(5,5,8), 0.01, 0.0, 3, 2, 0.1, 0.05) # default γ = 0.01

julia> result = [30.0 5.0 5.0 5.0 -20;
                 25.0 0.0 0.0 0.0 -25.0;
                 25.0 0.0 0.0 0.0 -25.0;
                 25.0 0.0 0.0 0.0 -25.0;
                 20.0 -5.0 -5.0 -5.0 -30.0];

julia> for i in eachindex(result)
           @test result[i] .≈ -100 .* pressure[i] atol=1e-12
       end
```

# References

- [Peschka et al.](https://www.pnas.org/content/116/19/9275)
- [Craster and Matar](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131)
- [Derjaguin and Churaev](https://www.sciencedirect.com/science/article/abs/pii/0021979778900565)

"""
function filmpressure!(output, f, dgrad, γ, θ, n, m, hmin, hcrit)
    hip, hjp, him, hjm, hipjp, himjp, himjm, hipjm = viewneighbors(dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, f, (1, 0))
    circshift!(hjp, f, (0, 1))
    circshift!(him, f, (-1, 0))
    circshift!(hjm, f, (0, -1))
    # Diagonal elements  
    circshift!(hipjp, f, (1, 1))
    circshift!(himjp, f, (-1, 1))
    circshift!(himjm, f, (-1, -1))
    circshift!(hipjm, f, (1, -1))
    #= 
    Disjoining pressure part:
    1. Constant part due to angle, n, m, hmin
    2. Part due to the powerlaw
    =#
    if n == 9 && m == 3
        output .=
            -γ .* (
                (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .*
                (Swalbe.fast_93.(hmin ./ (f .+ hcrit)))
            )
    elseif n == 3 && m == 2
        output .=
            -γ .* (
                (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .*
                (Swalbe.fast_32.(hmin ./ (f .+ hcrit)))
            )
    else
        throw(
            DomainError(
                (n, m),
                "These exponents have not been used so far please take a look at `pressure.jl` and open an issue if there are questions",
            ),
        )
    end
    output .-=
        γ .* (
            2 / 3 .* (hjp .+ hip .+ him .+ hjm) .+
            1 / 6 .* (hipjp .+ himjp .+ himjm .+ hipjm) .- 10 / 3 .* f
        )
    return nothing
end
# Film pressure with the state struct
function filmpressure!(
    state::LBM_state_2D,
    sys::SysConst;
    θ = sys.param.θ,
    γ = sys.param.γ,
    n = sys.param.n,
    m = sys.param.m,
    hmin = sys.param.hmin,
    hcrit = sys.param.hcrit,
)
    hip, hjp, him, hjm, hipjp, himjp, himjm, hipjm = viewneighbors(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height, (1, 0))
    circshift!(hjp, state.height, (0, 1))
    circshift!(him, state.height, (-1, 0))
    circshift!(hjm, state.height, (0, -1))
    # Diagonal elements  
    circshift!(hipjp, state.height, (1, 1))
    circshift!(himjp, state.height, (-1, 1))
    circshift!(himjm, state.height, (-1, -1))
    circshift!(hipjm, state.height, (1, -1))
    # First the contact angle parameter part
    state.pressure .=
        -γ .* (
            (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .* (
                Swalbe.power_broad.(hmin ./ (state.height .+ hcrit), n) .-
                Swalbe.power_broad.(hmin ./ (state.height .+ hcrit), m)
            )
        )
    # Now the gradient
    state.pressure .-=
        γ .* (
            2 / 3 .* (hjp .+ hip .+ him .+ hjm) .+
            1 / 6 .* (hipjp .+ himjp .+ himjm .+ hipjm) .- 10 / 3 .* state.height
        )
    return nothing
end

function filmpressure!(
    state::Expanded_2D,
    sys::SysConst;
    θ = sys.param.θ,
    γ = sys.param.γ,
    n = sys.param.n,
    m = sys.param.m,
    hmin = sys.param.hmin,
    hcrit = sys.param.hcrit,
)
    hip, hjp, him, hjm, hipjp, himjp, himjm, hipjm = viewneighbors(state.basestate.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.basestate.height, (1, 0))
    circshift!(hjp, state.basestate.height, (0, 1))
    circshift!(him, state.basestate.height, (-1, 0))
    circshift!(hjm, state.basestate.height, (0, -1))
    # Diagonal elements  
    circshift!(hipjp, state.basestate.height, (1, 1))
    circshift!(himjp, state.basestate.height, (-1, 1))
    circshift!(himjm, state.basestate.height, (-1, -1))
    circshift!(hipjm, state.basestate.height, (1, -1))
    # First the contact angle parameter part
    state.basestate.pressure .=
        -γ .* (
            (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .* (
                Swalbe.power_broad.(hmin ./ (state.basestate.height .+ hcrit), n) .-
                Swalbe.power_broad.(hmin ./ (state.basestate.height .+ hcrit), m)
            )
        )
    # Now the gradient
    state.basestate.pressure .-=
        γ .* (
            2 / 3 .* (hjp .+ hip .+ him .+ hjm) .+
            1 / 6 .* (hipjp .+ himjp .+ himjm .+ hipjm) .- 10 / 3 .* state.basestate.height
        )
    return nothing
end

# One dim implementation
function filmpressure!(output::Vector, f, dgrad, γ, θ, n, m, hmin, hcrit)
    hip, him = viewneighbors_1D(dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, f, 1)
    circshift!(him, f, -1)
    #= 
    Disjoining pressure part:
    1. Constant part due to angle, n, m, hmin
    2. Part due to the powerlaw
    =#
    if n == 9 && m == 3
        output .=
            -γ .* (
                (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .*
                (Swalbe.fast_93.(hmin ./ (f .+ hcrit)))
            )
    elseif n == 3 && m == 2
        output .=
            -γ .* (
                (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .*
                (Swalbe.fast_32.(hmin ./ (f .+ hcrit)))
            )
    else
        throw(
            DomainError(
                (n, m),
                "These exponents have not been used so far please take a look at `pressure.jl` and open an issue if there are questions",
            ),
        )
    end
    output .-= γ .* (hip .- 2 .* f .+ him)
    return nothing
end


function filmpressure!(
    state::LBM_state_1D,
    sys::Consts_1D;
    θ = sys.param.θ,
    n = sys.param.n,
    m = sys.param.m,
    hmin = sys.param.hmin,
    hcrit = sys.param.hcrit,
    γ = sys.param.γ,
)
    hip, him = viewneighbors_1D(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height, 1)
    circshift!(him, state.height, -1)

    @. state.pressure .=
        -γ * (
            (1 - cospi(θ)) * (n - 1) * (m - 1) / ((n - m) * hmin) * (
                Swalbe.power_broad(hmin / (state.height + hcrit), n) -
                Swalbe.power_broad(hmin / (state.height + hcrit), m)
            )
        )

    @. state.pressure .-= γ * (hip - 2 * state.height + him)
    return nothing
end

function filmpressure!(
    state::Expanded_1D,
    sys::Consts_1D;
    θ = sys.param.θ,
    n = sys.param.n,
    m = sys.param.m,
    hmin = sys.param.hmin,
    hcrit = sys.param.hcrit,
    γ = sys.param.γ,
)
    hip, him = viewneighbors_1D(state.basestate.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.basestate.height, 1)
    circshift!(him, state.basestate.height, -1)

    @. state.basestate.pressure .=
        -γ * (
            (1 - cospi(θ)) * (n - 1) * (m - 1) / ((n - m) * hmin) * (
                Swalbe.power_broad(hmin / (state.basestate.height + hcrit), n) -
                Swalbe.power_broad(hmin / (state.basestate.height + hcrit), m)
            )
        )

    @. state.basestate.pressure .-= γ * (hip - 2 * state.basestate.height + him)
    return nothing
end

function filmpressure!(
    state::State_gamma_1D,
    sys::Consts_1D;
    θ = sys.param.θ,
    n = sys.param.n,
    m = sys.param.m,
    hmin = sys.param.hmin,
    hcrit = sys.param.hcrit,
    γ = sys.param.γ,
)
    hip, him = viewneighbors_1D(state.basestate.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.basestate.height, 1)
    circshift!(him, state.basestate.height, -1)

    @. state.basestate.pressure .=
        -γ * (
            (1 - cospi(θ)) * (n - 1) * (m - 1) / ((n - m) * hmin) * (
                Swalbe.power_broad(hmin / (state.basestate.height + hcrit), n) -
                Swalbe.power_broad(hmin / (state.basestate.height + hcrit), m)
            )
        )
    # Should be fine as long as τ = 1
    ft0, ft1, ft2 = viewdists_1D(state.basestate.ftemp)
    # Save pressure contributions so one can evalute their overall contribution
    ft1 .= state.basestate.pressure
    ft2 .= -γ .* (hip .- 2 .* state.basestate.height .+ him)

    @. state.basestate.pressure .-= γ * (hip - 2 * state.basestate.height + him)
    return nothing
end

# Paolo active matter model
function filmpressure!(output::Vector, f, dgrad, rho, γ, θ, n, m, hmin, hcrit; Gamma = 0.0)
    hip, him = viewneighbors_1D(dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, f, 1)
    circshift!(him, f, -1)
    #= 
    Disjoining pressure part:
    1. Constant part due to angle, n, m, hmin
    2. Part due to the powerlaw
    =#
    output .=
        -(γ .+ Gamma .* rho) .* (
            (1 .- cospi.(θ)) .* (n - 1) .* (m - 1) ./ ((n - m) * hmin) .* (
                Swalbe.power_broad.(hmin ./ (f .+ hcrit), n) .-
                Swalbe.power_broad.(hmin ./ (f .+ hcrit), m)
            )
        )

    output .-= (γ .+ Gamma .* rho) .* (hip .- 2 .* f .+ him)
    return nothing
end


"""
	function filmpressure!(state::MultiLayer_2D, sys::SysConstMultiLayer)

Calculates the film pressure of a multilayer system. For 3D simulations we currently only have two layers implemented.

# Arguments

- `state`: the simulation fields
- `sys`: the system parmaeters

# Math

## Two layers

```math
            p_1-p_2 =&  \\varphi_{01}'(z_1)-\\varphi_{11}'(z_2-z_1)\\left(1+\\frac{|\\nabla z_1|^2}{2}\\right)-  (\\gamma_{12}+\\varphi_{11}(z_2-z_1)) \\nabla^2 z_1
```
```math
            p_2 - p_{atm} =& \\varphi_{11}'(z_2-z_1)\\left(1+\\frac{|\\nabla z_1|^2}{2}\\right)+\\varphi'_{02}(z_2)  - \\nabla\\cdot \\gamma_{2v} \\nabla z_2
```

## Three layers

```math
        p_1-p_2=&- \\left(\\gamma_{12}+\\varphi_{11}(h_2)+\\varphi_{12}(h_2+h_3)\\right) \\nabla^2 z_1 - \\left((\\varphi_{11}'(h_2)+\\varphi_{12}'(h_2+h_3)\\right)\\left(1+\\frac{(\\nabla z_1)^2}{2} \\right) + \\varphi'_{01}(h_1)\\left(1  + \\frac{(\\nabla z_0)^2}{2}\\right)
```
```math
       p_2-p_3=&-(\\gamma_{23}+\\varphi_{21}(h_3)) \\nabla^2 z_2- \\varphi'_{21}(h_3)\\left(1+ \\frac{(\\nabla z_2)^2}{2}\\right)+ \\varphi'_{02}(z_2)\\left(1+ \\frac{(\\nabla z_0)^2}{2}\\right)+\\varphi'_{11}(h_2)\\left(1+\\frac{(\\nabla z_1)^2}{2}\\right)
```
```math
         p_3-p_{atm}=&-\\gamma_{3v}\\nabla^2 z_3 + \\varphi'_{03}(z_3)\\left(1 + \\frac{(\\nabla z_0)^2}{2}\\right)  + \\varphi'_{12}(h_2+h_3)\\left(1  + \\frac{(\\nabla z_1)^2}{2}\\right)+ \\varphi'_{21}(h_3)\\left(1 + \\frac{(\\nabla z_2)^2}{2}\\right).
```

# Reference

- [Richter et al.](https://arxiv.org/abs/2409.16659)
"""
function filmpressure!(state::MultiLayer_2D, sys::SysConstMultiLayer)
    #disjoining pressure, Remark that \phi'(h_1+h_2) comes with a factor of 2 in front of hmin and hcrit
    ∇f!(state.grad_hx, state.grad_hy, state.height[:,:,1], state.dgrad[:,:,1,:], 1)
    state.grad_h_sq .= state.grad_hx .* state.grad_hx .+ state.grad_hy .* state.grad_hy
      #Laplace pressure
    hip, hjp, him, hjm, hipjp, himjp, himjm, hipjm = viewneighborsMultiLayer(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height, (1,0))
    circshift!(hjp, state.height, (0,1))
    circshift!(him, state.height, (-1,0))
    circshift!(hjm, state.height, (0,-1))
    # Diagonal elements
    circshift!(hipjp, state.height, (1,1))
    circshift!(himjp, state.height, (-1,1))
    circshift!(himjm, state.height, (-1,-1))
    circshift!(hipjm, state.height, (1,-1))
	if sys.n==9 && sys.m==3
       # if sys.layers == 2
       # store the one field we need twice
       state.hi[:,:,1].=(sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)   .* Swalbe.fast_93.(sys.hmin./(state.height[:,:,2] .+ sys.hcrit)).* (1 .+ state.grad_h_sq ./2)
        # disjoining pressure
        state.pressure[:,:,1] .=  (
                (sys.gamma[1,3]-sys.gamma[1,2]-sys.gamma[2,3])                  *(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)   .* Swalbe.fast_93.(sys.hmin./(state.height[:,:,1] .+ sys.hcrit))
                .- state.hi[:,:,1]
            )
            state.pressure[:,:,2] .=  (
                state.hi[:,:,1]
                .+ (sys.gamma[2,3]+sys.gamma[1,4]-sys.gamma[2,4]-sys.gamma[1,3])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)   .* Swalbe.fast_93.(2*sys.hmin./(state.height[:,:,1] .+ state.height[:,:,2] .+ 2*sys.hcrit))
            )
            #Laplace pressure
            state.pressure[:,:,1] .-= (sys.gamma[2,3] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4]) .* 1/6 .* Swalbe.wetting_potential_93.(sys.hmin./(state.height[:,:,2] .+ sys.hcrit))) .* (
                     2/3 .* (hjp[:,:,1] .+ hip[:,:,1] .+ him[:,:,1] .+ hjm[:,:,1])
                     .+ 1/6 .* (hipjp[:,:,1] .+ himjp[:,:,1] .+ himjm[:,:,1] .+ hipjm[:,:,1])
                     .- 10/3 .* state.height[:,:,1]
            )
            state.pressure[:,:,2] .-= (
                (
                    (
                        2/3 .* (hjp[:,:,1] .+ hip[:,:,1] .+ him[:,:,1] .+ hjm[:,:,1])
                       .+ 1/6 .* (hipjp[:,:,1] .+ himjp[:,:,1] .+ himjm[:,:,1] .+ hipjm[:,:,1])
                       .- 10/3 .* state.height[:,:,1]
                    ).+(
                        2/3 .* (hjp[:,:,2] .+ hip[:,:,2] .+ him[:,:,2] .+ hjm[:,:,2])
                        .+ 1/6 .* (hipjp[:,:,2] .+ himjp[:,:,2] .+ himjm[:,:,2] .+ hipjm[:,:,2])
                        .- 10/3 .* state.height[:,:,2]
                    )
                ).*  sys.gamma[3,4]
             )
        # end
	elseif sys.n==3 && sys.m==2
      # if sys.layers == 2
      #Store what we will need twice
      state.hi[:,:,1] .= (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)   .* Swalbe.fast_32.(sys.hmin  ./(state.height[:,:,2] .+ sys.hcrit)) .* (1 .+ state.grad_h_sq ./2)
        # disjoining pressure
        state.pressure[:,:,1] .=  (
            (sys.gamma[1,3]-sys.gamma[1,2]-sys.gamma[2,3])                  *(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)   .* Swalbe.fast_32.(sys.hmin./(state.height[:,:,1] .+ sys.hcrit))
            .- state.hi[:,:,1]
        )
        state.pressure[:,:,2] .=  (
            state.hi[:,:,1]
            .+ (sys.gamma[2,3]+sys.gamma[1,4]-sys.gamma[2,4]-sys.gamma[1,3])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)   .* Swalbe.fast_32.(2*sys.hmin  ./(state.height[:,:,1] .+ state.height[:,:,2] + 2*sys.hcrit))
        )
        #Laplace pressure
        state.pressure[:,:,1] .-= (
                (
                     2/3 .* (hjp[:,:,1] .+ hip[:,:,1] .+ him[:,:,1] .+ hjm[:,:,1])
                     .+ 1/6 .* (hipjp[:,:,1] .+ himjp[:,:,1] .+ himjm[:,:,1] .+ hipjm[:,:,1])
                     .- 10/3 .* state.height[:,:,1]
                ).* (sys.gamma[2,3] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4]) .* 1/6 .* Swalbe.wetting_potential_32.(sys.hmin./(state.height[:,:,2] .+ sys.hcrit)))
            )
            state.pressure[:,:,2] .-= (
                    (
                        2/3 .* (hjp[:,:,1] .+ hip[:,:,1] .+ him[:,:,1] .+ hjm[:,:,1])
                       .+ 1/6 .* (hipjp[:,:,1] .+ himjp[:,:,1] .+ himjm[:,:,1] .+ hipjm[:,:,1])
                       .- 10/3 .* state.height[:,:,1]
                    ).+(
                        2/3 .* (hjp[:,:,2] .+ hip[:,:,2] .+ him[:,:,2] .+ hjm[:,:,2])
                        .+ 1/6 .* (hipjp[:,:,2] .+ himjp[:,:,2] .+ himjm[:,:,2] .+ hipjm[:,:,2])
                        .- 10/3 .* state.height[:,:,2]
                    )
                ).*  sys.gamma[3,4]
    # end
	else
        	throw(DomainError((sys.n,sys.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
    	end
    return nothing
end




"""
    surface_tension!(state::StateMiscible_1D, sys::SysConstMiscible_1D)

Calculates and updates the surface tensions for a miscible two-phase system in 1D based on the weighted average of individual surface tensions of each phase.

# Parameters

- `state::StateMiscible_1D`:
   An object that holds the current state of the system.
   - `state.height[:,1]` and `state.height[:,2]` represent the heights of the two phases in each cell.
   - `state.gamma[:,i]` (i = 1, 2, 3) stores the surface tension values, i=1 for the liquid-gas surface tension, i=2 for the solid-liquid surface tension, i3 for the solid-vapour surface tensions
- `sys::SysConstMiscible_1D`:
   An object that contains system constants.
   - `sys.gamma[i,j]` represents the surface tension for the jth liquid, between solid-vapour (i=1), solid-liquid (i=2), liquid-vapour (i=3)

# Calculations

For each spatial component (denoted by indices `i = 1, 2, 3`):
- Computes the effective surface tension as the height-weighted average:
  \\[
  \\gamma_{i} = \\frac{h_1 \\gamma_{i,1} + h_2 \\gamma_{i,2}}{h_1 + h_2},
  \\]
  where:
  - \\( h_1, h_2 \\) are the heights of the two phases,
  - \\( \\gamma_{i,1}, \\gamma_{i,2} \\) are the surface tension constants and i in (sv, sl, lv).

# Notes
- The function updates the `gamma` values in-place for efficiency.
- Ensure that `state.height` does not contain zeros to avoid division by zero.

# Example Usage
```julia
surface_tension!(state, sys)
```
After execution, `state.gamma` will hold the updated surface tension values for each spatial component.
"""
function surface_tension!(state::StateMiscible_1D, sys::SysConstMiscible_1D)
    state.gamma[:,1] .= (state.height[:,1] .* sys.gamma[1,1] .+ state.height[:,2] .* sys.gamma[1,2])./(state.height[:,1] .+ state.height[:,2])
    state.gamma[:,2] .= (state.height[:,1] .* sys.gamma[2,1] .+ state.height[:,2] .* sys.gamma[2,2])./(state.height[:,1] .+ state.height[:,2])
    state.gamma[:,3] .= (state.height[:,1] .* sys.gamma[3,1] .+ state.height[:,2] .* sys.gamma[3,2])./(state.height[:,1] .+ state.height[:,2])
end


#This moving averadge stabilizes the simulation of droplet coalescence against a surface tension gradient but destabilizes the solvation of a droplet with high surface tension.

"""
    surface_tension_smooth!(state::StateMiscible_1D, sys::SysConstMiscible_1D; center_weight=0.5)

Calculates and updates the surface tensions for a miscible two-phase system in 1D, then applies a smoothing operation to stabilize simulations against surface tension gradients.

# Parameters
- `state::StateMiscible_1D`:
   An object that holds the current state of the system.
   - `state.height[:,1]` and `state.height[:,2]` represent the heights of the two phases in each cell.
   - `state.gamma[:,i]` (i = 1, 2, 3) stores the surface tension values:
     - \\( i=1 \\): liquid-gas surface tension,
     - \\( i=2 \\): solid-liquid surface tension,
     - \\( i=3 \\): solid-vapor surface tension.
- `sys::SysConstMiscible_1D`:
   An object that contains system constants.
   - `sys.gamma[i,j]` represents the surface tension for the \\( j \\)-th liquid, between:
     - solid-vapor (\\( i=1 \\)),
     - solid-liquid (\\( i=2 \\)),
     - liquid-vapor (\\( i=3 \\)).
- `center_weight=0.5` (optional):
   The weight assigned to the original surface tension value when applying the smoothing operation. Defaults to 0.5.

# Calculations
1. **Initial Surface Tension Calculation**:
   For each spatial component \\( i \\) (1, 2, 3), the effective surface tension is computed as:
   \\[
   \\gamma_{i} = \\frac{h_1 \\gamma_{i,1} + h_2 \\gamma_{i,2}}{h_1 + h_2},
   \\]
   where:
   - \\( h_1, h_2 \\) are the heights of the two phases,
   - \\( \\gamma_{i,1}, \\gamma_{i,2} \\) are the surface tension constants.

2. **Smoothing**:
   - For each spatial component \\( \\gamma_{i} \\), the smoothed value is computed using neighbor contributions:
     \\[
     \\gamma_{i,\\text{smoothed}} = w \\gamma_{i,\\text{original}} + (1 - w)(\\gamma_{i,\\text{left}} + \\gamma_{i,\\text{right}}),
     \\]
     where \\( w \\) is the `center_weight`, and `left`/`right` neighbors are determined via `circshift!`.

# Purpose
This smoothing process enhances numerical stability in simulations involving droplet coalescence under surface tension gradients. However, it may destabilize the solvation of droplets with high surface tensions.

# Notes
- The function updates the `gamma` values in-place for efficiency.
- Ensure that `state.height` does not contain zeros to avoid division by zero.
- Neighboring operations use `viewneighbors_1D(state.dgrad[:,:,1])`, so ensure compatibility of `state` with this utility.

# Example Usage
```julia
surface_tension_smooth!(state, sys; center_weight=0.7)
```
This applies smoothing with a higher weighting on the original value.
"""

function surface_tension_smooth!(state::StateMiscible_1D, sys::SysConstMiscible_1D; center_weight=0.5)
    state.gamma[:,1] .= (state.height[:,1] .* sys.gamma[1,1] .+ state.height[:,2] .* sys.gamma[1,2])./(state.height[:,1] .+ state.height[:,2])
    state.gamma[:,2] .= (state.height[:,1] .* sys.gamma[2,1] .+ state.height[:,2] .* sys.gamma[2,2])./(state.height[:,1] .+ state.height[:,2])
    state.gamma[:,3] .= (state.height[:,1] .* sys.gamma[3,1] .+ state.height[:,2] .* sys.gamma[3,2])./(state.height[:,1] .+ state.height[:,2])
    gip, gim = viewneighbors_1D(state.dgrad[:,:,1])
    circshift!(gip, state.gamma[:,1], 1)
    circshift!(gim, state.gamma[:,1],-1)
    state.gamma[:,1] .= center_weight .* state.gamma[:,1] .+ (1-center_weight) .* (gip .+ gim)
    circshift!(gip, state.gamma[:,2], 1)
    circshift!(gim, state.gamma[:,2],-1)
    state.gamma[:,2] .= center_weight .* state.gamma[:,2] .+ (1-center_weight) .* (gip .+ gim)
    circshift!(gip, state.gamma[:,3], 1)
    circshift!(gim, state.gamma[:,3],-1)
    state.gamma[:,3] .= center_weight .* state.gamma[:,3] .+ (1-center_weight) .* (gip .+ gim)
end


function filmpressure!(state::StateMiscible_1D, sys::SysConstMiscible_1D)
    hip, him = viewneighbors_1D(state.dgrad[:,:,1])
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height[:,1].+state.height[:,2], 1)
    circshift!(him, state.height[:,1].+state.height[:,2], -1)
	if sys.n==9 && sys.m==3
    		state.pressure[:,1] .= (-state.gamma[:,1] .- state.gamma[:,2] .+ state.gamma[:,3]) .* (sys.n - 1) .* (sys.m - 1) ./ ((sys.n - sys.m) * sys.hmin) .* Swalbe.fast_93.(sys.hmin./(state.height[:,1] .+ state.height[:,2] .+ 2 .* sys.hcrit))
            state.pressure[:,2] .=state.pressure[:,1]
            # For extra stability
	    state.pressure .+= 0.1.* (-state.gamma[:,1] .- state.gamma[:,2] .+ state.gamma[:,3]) .* (sys.n - 1) .* (sys.m - 1) ./ ((sys.n - sys.m) * sys.hmin / 2) .* Swalbe.power_9.((sys.hmin/2)./(state.height .+ sys.hcrit))
	elseif sys.n==3 && sys.m==2
    		state.pressure .= (-state.gamma[:,1] .- state.gamma[:,2] .+ state.gamma[:,3]) .* (sys.n - 1) .* (sys.m - 1) ./ ((sys.n - sys.m) * sys.hmin) .* Swalbe.fast_32.(sys.hmin./(state.height[:,1] .+ state.height[:,2] .+ sys.hcrit))
	else
        	throw(DomainError((sys.n,sys.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
    	end
    state.pressure .-= state.gamma[:,1] .* (hip .- 2 .* (state.height[:,1] .+ state.height[:,2]) .+ him)
    return nothing
end


function filmpressure_fast!(state::Active_1D, sys::SysConstActive_1D)
    # All moved to the function surface_tension in active.jl
    # state.γ .= sys.γ_0 .- sys.Γ .* state.rho
    state.k .= ((-state.γ  .+ cospi.(sys.θ_0)*sys.γ_0 ) .* (sys.n - 1) .* (sys.m - 1) ./ ((sys.n - sys.m) * sys.hmin))
    #Π(h)
    if sys.n==9 && sys.m==3
        state.pressure .= (state.k .* Swalbe.fast_93.(sys.hmin ./ (state.height .+ sys.hcrit)))
    elseif sys.n==3 && sys.m==2
        state.pressure .= (state.k .* Swalbe.fast_32.(sys.hmin ./ (state.height .+ sys.hcrit)))
    else
        throw(DomainError((sys.n,sys.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
    end
    hip, him = viewneighbors_1D(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height, 1)
    circshift!(him, state.height, -1)
    state.pressure .-= state.γ .* (hip .- 2 .* state.height .+ him)
    return nothing
end


"""
	function filmpressure!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)

Calculates the film pressure of a multilayer system in pseudo 2D. Currently we have implemented two and three layers

# Arguments

- `state`: the simulation fields
- `sys`: the system parmaeters

# Math

## Two layers

```math
            p_1-p_2 =&  \\varphi_{01}'(z_1)-\\varphi_{11}'(z_2-z_1)\\left(1+\\frac{|\\nabla z_1|^2}{2}\\right)-  (\\gamma_{12}+\\varphi_{11}(z_2-z_1)) \\nabla^2 z_1
```
```math
            p_2 - p_{atm} =& \\varphi_{11}'(z_2-z_1)\\left(1+\\frac{|\\nabla z_1|^2}{2}\\right)+\\varphi'_{02}(z_2)  - \\nabla\\cdot \\gamma_{2v} \\nabla z_2
```

## Three layers

```math
        p_1-p_2=&- \\left(\\gamma_{12}+\\varphi_{11}(h_2)+\\varphi_{12}(h_2+h_3)\\right) \\nabla^2 z_1 - \\left((\\varphi_{11}'(h_2)+\\varphi_{12}'(h_2+h_3)\\right)\\left(1+\\frac{(\\nabla z_1)^2}{2} \\right) + \\varphi'_{01}(h_1)\\left(1  + \\frac{(\\nabla z_0)^2}{2}\\right)
```
```math
       p_2-p_3=&-(\\gamma_{23}+\\varphi_{21}(h_3)) \\nabla^2 z_2- \\varphi'_{21}(h_3)\\left(1+ \\frac{(\\nabla z_2)^2}{2}\\right)+ \\varphi'_{02}(z_2)\\left(1+ \\frac{(\\nabla z_0)^2}{2}\\right)+\\varphi'_{11}(h_2)\\left(1+\\frac{(\\nabla z_1)^2}{2}\\right)
```
```math
         p_3-p_{atm}=&-\\gamma_{3v}\\nabla^2 z_3 + \\varphi'_{03}(z_3)\\left(1 + \\frac{(\\nabla z_0)^2}{2}\\right)  + \\varphi'_{12}(h_2+h_3)\\left(1  + \\frac{(\\nabla z_1)^2}{2}\\right)+ \\varphi'_{21}(h_3)\\left(1 + \\frac{(\\nabla z_2)^2}{2}\\right).
```

# Reference

- [Richter et al.](https://arxiv.org/abs/2409.16659)
"""
function filmpressure!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)
    # We will need some derivatives of h
    # Maybe it is possible to do this more elegant without the for loop
    if sys.layers==2
        ∇f!(state.grad_h[:,1], state.height[:,1], state.dgrad[:,1,:])
    elseif sys.layers==3
        ∇f!(state.grad_h[:,1], state.height[:,1], state.dgrad[:,1,:])
        ∇f!(state.grad_h[:,2], state.height[:,1] .+ state.height[:,2], state.dgrad[:,1,:])
    end
    # And we will need some shifted values for the Laplacian later on
    hip, him = viewneighborsMultiLayer_1D(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height, 1)
    circshift!(him, state.height, -1)

    #### Two layer system
    if sys.layers==2
        if sys.n==9 && sys.m==3
             #Laplace pressure
            state.pressure[:,1] .= - (hip[:,1] .- 2 .* state.height[:,1] .+ him[:,1])                                                .* (sys.gamma[2,3] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4]) .* 1/6 .* Swalbe.wetting_potential_93.(sys.hmin./(state.height[:,2] .+ sys.hcrit)))
            state.pressure[:,2] .= - (hip[:,1] .+ hip[:,2] .- 2 .* (state.height[:,1] .+ state.height[:,2]) .+ him[:,1] .+ him[:,2]) .*  sys.gamma[3,4]
            # store the one field we need twice
            state.hi[:,1].=(sys.gamma[2,3]+sys.gamma[1,4]-sys.gamma[2,4]-sys.gamma[1,3])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*2*sys.hmin)   .* Swalbe.fast_93.(2*sys.hmin./(state.height[:,1] .+ state.height[:,2] .+ 2*sys.hcrit))
            # disjoining pressure
            state.pressure[:,1] .+= state.pressure[:,2] .+ state.hi[:,1] .+ (sys.gamma[1,3]-sys.gamma[1,2]-sys.gamma[2,3])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin) .* Swalbe.fast_93.(sys.hmin./(state.height[:,1] .+ sys.hcrit))
            state.pressure[:,2] .+=                        state.hi[:,1] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin) .* Swalbe.fast_93.(sys.hmin./(state.height[:,2] .+ sys.hcrit)).* (1 .+ state.grad_h[:,1] .* state.grad_h[:,1] ./2)
        elseif sys.n==3 && sys.m==2
           #Laplace pressure
           state.pressure[:,1] .= - (hip[:,1] .- 2 .* state.height[:,1] .+ him[:,1])                                                .* (sys.gamma[2,3] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4]) .* Swalbe.wetting_potential_32.(sys.hmin./(state.height[:,2] .+ sys.hcrit)))
           state.pressure[:,2] .= - (hip[:,1] .+ hip[:,2] .- 2 .* (state.height[:,1] .+ state.height[:,2]) .+ him[:,1] .+ him[:,2]) .*  sys.gamma[3,4]
       # store the one field we need twice
       state.hi[:,1].=(sys.gamma[2,3]+sys.gamma[1,4]-sys.gamma[2,4]-sys.gamma[1,3])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*2*sys.hmin)   .* Swalbe.fast_32.(2*sys.hmin./(state.height[:,1] .+ state.height[:,2] .+ 2*sys.hcrit))
           # disjoining pressure
           state.pressure[:,1] .+= state.pressure[:,2] .+ state.hi[:,1] .+ (sys.gamma[1,3]-sys.gamma[1,2]-sys.gamma[2,3])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin) .* Swalbe.fast_32.(  sys.hmin./(state.height[:,1] .+ sys.hcrit))
           state.pressure[:,2] .+=                        state.hi[:,1] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4])*(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin) .* Swalbe.fast_32.(sys.hmin./(state.height[:,2] .+ sys.hcrit)).* (1 .+ state.grad_h[:,1] .* state.grad_h[:,1] ./2)
        else
                throw(DomainError((sys.n,sys.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
        end

    ##### Three layer system
    elseif sys.layers==3
            if sys.n==3 && sys.m==2
                # Laplace pressure
                state.pressure[:,1] .= - (hip[:,1] .- 2 .* state.height[:,1] .+ him[:,1])                                                                                             .* (sys.gamma[2,3] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4]) .* Swalbe.wetting_potential_32.(sys.hmin./(state.height[:,2] .+ sys.hcrit)) .+ (sys.gamma[1,5] + sys.gamma[3,4] - sys.gamma[2,4] - sys.gamma[3,5])  .* Swalbe.wetting_potential_32.(2*sys.hmin ./ (state.height[:,2] .+ state.height[:,3] .+ 2*sys.hcrit)))
                state.pressure[:,2] .= - (hip[:,1] .+ hip[:,2] .- 2 .* (state.height[:,1] .+ state.height[:,2]) .+ him[:,1] .+ him[:,2])                                              .* (sys.gamma[3,4] .+ (sys.gamma[3,5]-sys.gamma[3,4]-sys.gamma[4,5]) .* Swalbe.wetting_potential_32.(sys.hmin./(state.height[:,3] .+ sys.hcrit)))
                state.pressure[:,3] .= - (hip[:,1] .+ hip[:,2]  .+ hip[:,3].- 2 .* (state.height[:,1] .+ state.height[:,2] .+ state.height[:,3]) .+ him[:,1] .+ him[:,2] .+ him[:,3]) .*  sys.gamma[4,5]
            # We will use this quite a lot so it cleans up a bit having it here
            prefac=(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)
            # Store things we will need multiple times
            # When going to more then 3 layers one should really think about adding an extra field to store those especially because then we could use state.hi to store the actual interface positions
                # \phi_123'(h_2)
                state.hi[:,1] .= (sys.gamma[1,5] + sys.gamma[2,4] - sys.gamma[1,4] - sys.gamma[2,5]) * 1/3 * prefac * Swalbe.fast_32.(3*sys.hmin ./ (state.height[:,1] .+ state.height[:,2] .+ state.height[:,3] .+ 3*sys.hcrit))
                # \phi_{23}'(h_2+h_3)(1+(\nabla z_1)^2/2)
                state.hi[:,2] .= (sys.gamma[2,5] + sys.gamma[3,4] - sys.gamma[2,4] -sys.gamma[3,5]) * 0.5 * prefac .* Swalbe.fast_32.( 2*sys.hmin ./ (state.height[:,2] .+ state.height[:,3] .+ 2* sys.hcrit)) .* (1 .+ state.grad_h[:,1] .* state.grad_h[:,1] ./ 2)
                # \phi_12'(h_3)(1+ (\nabla z_2)^2/2)
                state.hi[:,3] .=(sys.gamma[1,4] + sys.gamma[2,3] - sys.gamma[1,3] - sys.gamma[2,4]) * 0.5 * prefac * Swalbe.fast_32.(2*sys.hmin ./ (state.height[:,1] .+ state.height[:,2] .+ 2*sys.hcrit))
            # Disjoining pressure
                state.pressure[:,1] .+= state.pressure[:,2] .+ state.pressure[:,3] .+ state.hi[:,1]                  .+ state.hi[:,3] .+ (sys.gamma[1,3] - sys.gamma[1,2] - sys.gamma[2,3]) * prefac .* Swalbe.fast_32.(sys.hmin ./ (state.height[:,1] .+ sys.hcrit))
                state.pressure[:,2] .+= state.pressure[:,3]                        .+ state.hi[:,1] .+ state.hi[:,2] .+ state.hi[:,3] .+ (sys.gamma[2,4] - sys.gamma[2,3] - sys.gamma[3,4]) * prefac .* Swalbe.fast_32.(sys.hmin ./ (state.height[:,2] .+ sys.hcrit)) .* (1 .+ state.grad_h[:,1] .* state.grad_h[:,1] ./ 2)
                state.pressure[:,3] .+=                                               state.hi[:,1] .+ state.hi[:,2]                  .+ (sys.gamma[3,5] - sys.gamma[3,4] - sys.gamma[4,5]) * prefac .* Swalbe.fast_32.(sys.hmin ./ (state.height[:,3] .+ sys.hcrit)) .* (1 .+ state.grad_h[:,2] .* state.grad_h[:,2] ./ 2)
            elseif sys.n==9 && sys.m==3
                 # Laplace pressure
                 state.pressure[:,1] .= - (hip[:,1] .- 2 .* state.height[:,1] .+ him[:,1])                                                                                             .* (sys.gamma[2,3] .+ (sys.gamma[2,4]-sys.gamma[2,3]-sys.gamma[3,4]) .* 1/6 .* Swalbe.wetting_potential_93.(sys.hmin./(state.height[:,2] .+ sys.hcrit)) .+ (sys.gamma[1,5] + sys.gamma[3,4] - sys.gamma[2,4] - sys.gamma[3,5]) * 1/6 .* Swalbe.wetting_potential_93.(2*sys.hmin ./ (state.height[:,2] .+ state.height[:,3] .+ 2*sys.hcrit)))
                 state.pressure[:,2] .= - (hip[:,1] .+ hip[:,2] .- 2 .* (state.height[:,1] .+ state.height[:,2]) .+ him[:,1] .+ him[:,2])                                              .* (sys.gamma[3,4] .+ (sys.gamma[3,5]-sys.gamma[3,4]-sys.gamma[4,5]) .* 1/6 .* Swalbe.wetting_potential_93.(sys.hmin./(state.height[:,3] .+ sys.hcrit)))
                 state.pressure[:,3] .= - (hip[:,1] .+ hip[:,2]  .+ hip[:,3].- 2 .* (state.height[:,1] .+ state.height[:,2] .+ state.height[:,3]) .+ him[:,1] .+ him[:,2] .+ him[:,3]) .*  sys.gamma[4,5]
             # We will use this quite a lot so it cleans up a bit having it here
             prefac=(sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)
             # Store things we will need multiple times
             # When going to more then 3 layers one should really think about adding an extra field to store those especially because then we could use state.hi to store the actual interface positions
                 # \phi_123'(h_2)
                 state.hi[:,1] .= (sys.gamma[1,5] + sys.gamma[2,4] - sys.gamma[1,4] - sys.gamma[2,5]) * 1/3 * prefac *Swalbe.fast_93.(3*sys.hmin ./ (state.height[:,1] .+ state.height[:,2] .+ state.height[:,3] .+ 3*sys.hcrit))
                 # \phi_{23}'(h_2+h_3)(1+(\nabla z_1)^2/2)
                 state.hi[:,2] .= (sys.gamma[2,5] + sys.gamma[3,4] - sys.gamma[2,4] -sys.gamma[3,5]) * 0.5 * prefac .*Swalbe.fast_93.( 2*sys.hmin ./ (state.height[:,2] .+ state.height[:,3] .+ 2* sys.hcrit)) .* (1 .+ state.grad_h[:,1] .* state.grad_h[:,1] ./ 2)
                 # \phi_12'(h_3)(1+ (\nabla z_2)^2/2)
                 state.hi[:,3] .=(sys.gamma[1,4] + sys.gamma[2,3] - sys.gamma[1,3] - sys.gamma[2,4]) * 0.5 * prefac *Swalbe.fast_93.(2*sys.hmin ./ (state.height[:,1] .+ state.height[:,2] .+ 2*sys.hcrit))
             # Disjoining pressure
                 state.pressure[:,1] .+= state.pressure[:,2] .+ state.pressure[:,3] .+ state.hi[:,1]                  .+ state.hi[:,3] .+ (sys.gamma[1,3] - sys.gamma[1,2] - sys.gamma[2,3]) * prefac .*Swalbe.fast_93.(sys.hmin ./ (state.height[:,1] .+ sys.hcrit))
                 state.pressure[:,2] .+= state.pressure[:,3]                        .+ state.hi[:,1] .+ state.hi[:,2] .+ state.hi[:,3] .+ (sys.gamma[2,4] - sys.gamma[2,3] - sys.gamma[3,4]) * prefac .*Swalbe.fast_93.(sys.hmin ./ (state.height[:,2] .+ sys.hcrit)) .* (1 .+ state.grad_h[:,1] .* state.grad_h[:,1] ./ 2)
        else
                throw(DomainError((sys.n,sys.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
        end
    end
    return nothing
end


"""
    function extra_pressure!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)

Attempt to keep the precursor films up
"""
function extra_pressure!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)
   # And we will need some shifted values for the Laplacian later on
   hip, him = viewneighborsMultiLayer_1D(state.dgrad)
   # Straight elements j+1, i+1, i-1, j-1
   circshift!(hip, state.height, 1)
   circshift!(him, state.height, -1)
   state.pressure[:,1] .-=  sys.gamma[2,3] .* (hip[:,1] .- 2 .* state.height[:,1] .+ him[:,1]) .* sys.extra_pressure_fac .* Swalbe.power_6.(sys.hmin ./ (state.height[:,1] .+ sys.hcrit))
   state.pressure[:,2] .-=  sys.gamma[3,4] .* (hip[:,2] .- 2 .* state.height[:,2] .+ him[:,2]) .* sys.extra_pressure_fac .* Swalbe.power_6.(sys.hmin ./ (state.height[:,2] .+ sys.hcrit))
   state.pressure[:,3] .-=  sys.gamma[4,5] .* (hip[:,3] .- 2 .* state.height[:,3] .+ him[:,3]) .* sys.extra_pressure_fac .* Swalbe.power_6.(sys.hmin ./ (state.height[:,3] .+ sys.hcrit))
end

"""
	filmpressure_curved!(state::State_curved_1D, sys::SysConst_1D)

Computes the capillary pressure for a height field, taking into account the curvature of an underlying substrate.
This extended model incorporates an additional term to account for the influence of a curved substrate in
the capillary pressure computation, building on the base functionality of `filmpressure!` by introducing a
substrate curvature-dependent factor. Compare with [`filmpressure!`](@ref).

# Arguments
- `state::State_curved_1D`: A struct containing the state variables for the current 1D problem.
    - `state.height`: The current height field `h(x,t)`.
    - `state.substrate`: The substrate shape field.
    - `state.grad_substrate`: Gradient of the substrate, pre-computed to account for curved geometry.
    - `state.dgrad`: Discretization step for numerical differentiation.
    - `state.pressure`: Array to store the resulting capillary pressure.
- `sys::SysConst_1D`: A struct containing system constants for the computation.
    - `sys.γ`: Surface tension coefficient.
    - `sys.θ`: Equilibrium contact angle.
    - `sys.n::Int`: Larger exponent for the disjoining pressure `Π(h)`.
    - `sys.m::Int`: Smaller exponent for the disjoining pressure `Π(h)`.
    - `sys.hmin`: Parameter of `Π(h)`, satisfying `Π(hmin) = 0`.
    - `sys.hcrit`: Stabilizer for the case `h(x,t) ≪ hmin`.

# Mathematics
The capillary pressure `p_cap` in this model accounts for the effects of substrate curvature and the disjoining pressure `Π(h)`.
It is expressed as:

```math
p_{\\text{cap}} = -(1 + 0.5 \\nabla s^2) \\gamma (1 - \\cos(\\theta)) \\phi'(h) - \\gamma \\nabla^2 (h+b)
```

# References
- Tilman Richter, Thin film modelling in the lubrication framework, Thesis
- [Peschka et al.](https://www.pnas.org/content/116/19/9275)
- [Craster and Matar](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131)
- [Derjaguin and Churaev](https://www.sciencedirect.com/science/article/abs/pii/0021979778900565)
"""
function filmpressure_curved!(state::State_curved_1D, sys::SysConst_1D)
    Swalbe.∇f!(state.grad_substrate, state.substrate, state.dgrad)
    hip, him = viewneighbors_1D(state.dgrad)
    # Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height .+ state.substrate, 1)
    circshift!(him, state.height .+ state.substrate, -1)
        if sys.param.n==9 && sys.param.m==3
                state.pressure .= -(1 .+ 0.5 .* state.grad_substrate .* state.grad_substrate) .* sys.param.γ .* (1 .- cospi.(sys.param.θ)) .* (sys.param.n - 1) .* (sys.param.m - 1) ./ ((sys.param.n - sys.param.m) * sys.param.hmin) .* Swalbe.fast_93.(sys.param.hmin./(state.height .+ sys.param.hcrit))
        elseif sys.param.n==3 && sys.param.m==2
                state.pressure .= - (1 .+ 0.5 .* state.grad_substrate .* state.grad_substrate) .* sys.param.γ .* (1 .- cospi.(sys.param.θ)) .* (sys.param.n - 1) .* (sys.param.m - 1) ./ ((sys.param.n - sys.param.m) * sys.param.hmin) .* Swalbe.fast_32.(sys.param.hmin./(state.height .+ sys.param.hcrit))
        elseif sys.param.n==0 && sys.param.m==0
                state.pressure .= 0.0
        else
                throw(DomainError((sys.param.n,sys.param.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
        end
        state.pressure .-= sys.param.γ .* (hip .- 2 .* (state.height .+ state.substrate) .+ him)
    return nothing
end



"""
	filmpressure_curved_healing!(state::State_curved_1D, sys::SysConst_1D, S)
	
Computes the capillary pressure for a height field on a curved substrate, allowing for a negative spreading coefficient `S`.
This enables the simulation of the healing length of a liquid film, which affects how the film spreads and stabilizes
on the substrate. This function builds on `filmpressure_curved!` but incorporates the additional influence of `S` on
the disjoining pressure term.

# Arguments
- `state::State_curved_1D`: Struct holding the current state variables in 1D.
    - `state.height`: The height field `h(x,t)`.
    - `state.substrate`: Shape field for the underlying substrate.
    - `state.grad_substrate`: Gradient of the substrate, pre-computed for efficient curvature effects.
    - `state.dgrad`: Step size for numerical gradient computations.                                                                                             - `state.pressure`: Array where the computed capillary pressure will be stored.
    - `sys.θ`: Equilibrium contact angle.
    - `sys.n::Int`: Larger exponent for the disjoining pressure `Π(h)`.
    - `sys.m::Int`: Smaller exponent for the disjoining pressure `Π(h)`.
    - `sys.hmin`: Parameter of `Π(h)`, satisfying `Π(hmin) = 0`.
    - `sys.hcrit`: Stabilizer for cases where `h(x,t) ≪ hmin`.
- `S::Number`: Spreading coefficient ``\\gamma_{sv}-(\\gamma+\\gamma_{sl})``

# Mathematics
The capillary pressure `p_cap` on a curved substrate is modified to include the spreading coefficient `S`, affecting the
disjoining pressure term:
```math
p_{\\text{cap}} = S (1 + 0.5 \\nabla s^2) \\phi'(h) - \\gamma \\nabla^2 h
```

# References
- Tilman Richter, Thin film modelling in the lubrication approximation                                                                                      - [Peschka et al.](https://www.pnas.org/content/116/19/9275)
- [Craster and Matar](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131)
- [Derjaguin and Churaev](https://www.sciencedirect.com/science/article/abs/pii/0021979778900565)                                                           """
function filmpressure_curved_healing!(state::State_curved_1D, sys::SysConst_1D,S)
    Swalbe.∇f!(state.grad_substrate, state.substrate, state.dgrad)
    hip, him = viewneighbors_1D(state.dgrad)
	# Straight elements j+1, i+1, i-1, j-1
    circshift!(hip, state.height .+ state.substrate, 1)
    circshift!(him, state.height .+ state.substrate, -1)
        if sys.param.n==9 && sys.param.m==3
                state.pressure .= S .* (1 .+ 0.5 .* state.grad_substrate .* state.grad_substrate)* (sys.param.n - 1) .* (sys.param.m - 1) ./ ((sys.param.n - sys.param.m) * sys.param.hmin) .* Swalbe.fast_93.(sys.param.hmin./(state.height .+ sys.param.hcrit))
        elseif sys.param.n==3 && sys.param.m==2                                                                                                                                         state.pressure .= S .* (1 .+ 0.5 .* state.grad_substrate .* state.grad_substrate) .* (sys.param.n - 1) .* (sys.param.m - 1) ./ ((sys.param.n - sys.param.m) * sys.param.hmin) .* Swalbe.fast_32.(sys.param.hmin./(state.height .+ sys.param.hcrit))
        elseif sys.param.n==0 && sys.param.m==0
                state.pressure .= 0.0
        else
		throw(DomainError((sys.param.n,sys.param.m), "This disjoining pressure is not implemented, Options currently are (n,m)=(9,3) or (n,m)=(3,2). Use those or implement a new option."))
	end
	state.pressure .-= sys.param.γ .* (hip .- 2 .* (state.height .+ state.substrate) .+ him)
    return nothing
end


"""
	function rho_pressure!(state::Active_1D, sys::SysConstActive_1D)

We apply an artificial pressure on the catalyst density that fullfils two roles

1. Makes sure that ``\\rho>0`` via
``p_1=-\\left(\\frac{h^*}{\\rho + h^*}\\right)^n.`
Without ``p_1`` there is advection of more catalyst then actually is present
2. Makes sure that catalyst does not diffuse into the precurso layer that we consider as dry via
``p_2=\\left(\\frac{h^*}{h + h^*}\\right)^n.``
The pressure for ``\\rho`` reads
``p=p_1+p_2``
"""
function rho_pressure!(state::Active_1D, sys::SysConstActive_1D)
    state.rho_pressure .=  (
				  .- Swalbe.power_6.(sys.hcrit ./(state.rho .+ sys.rho_crit))
				  .+ Swalbe.power_6.(sys.hcrit ./(state.height .+ sys.hcrit2))
			)
end


"""
    rho_A_pressure!(state::Active_1D, sys::SysConstActive_1D)

Applies an artificial pressure to `ρ_A`, the **product** density, ensuring that `ρ_A > 0` and preventing unintended diffusion of the product into the precursor layer.

# Arguments
- `state::Active_1D`: Contains `rho_A` (product density), `height`, and `rho_A_pressure`.
- `sys::SysConstActive_1D`: System constants, including `hcrit`, `b_A`, and `hcrit2`, used to define pressure terms for the product.

# Pressure Model

The pressure on `\\rho_A` is calculated as follows:
1. Ensures `\\rho_A > 0`:
   ``` p_1 = -\\left(\\frac{h^*}{\\rho_A + b_A}\\right)^6 ```
2. Prevents diffusion into the dry precursor layer:
   ``` p_2 = \\left(\\frac{h^*}{h + h^*}\\right)^6 ```
The total pressure for `\\rho_A` is given by `p = p_1 + p_2`.

# Returns

- `nothing`: Updates `state.rho_A_pressure` with the calculated artificial pressure.
"""
function rho_A_pressure!(state::Active_1D, sys::SysConstActive_1D)
    state.rho_A_pressure .=  (
				  .- Swalbe.power_6.(sys.hcrit ./(state.rho_A .+ sys.b_A))
				  .+ Swalbe.power_6.(sys.hcrit ./(state.height .+ sys.hcrit2))
			)
end



"""
    rho_B_pressure!(state::Active_1D, sys::SysConstActive_1D)

Calculates an artificial pressure on `ρ_B`, the **reactant** density, to ensure `ρ_B > 0` and prevent unwanted diffusion of the reactant into the precursor layer.

# Arguments
- `state::Active_1D`: Contains `rho_B` (reactant density), `height`, and `rho_B_pressure`.
- `sys::SysConstActive_1D`: System constants, including parameters `hcrit`, `b_B`, and `hcrit2` for the pressure terms affecting the reactant.

# Pressure Model

The pressure on `\\rho_B` is defined by two terms:
1. Maintains `\\rho_B > 0` with:
   ``` p_1 = -\\left(\\frac{h^*}{\\rho_B + b_B}\\right)^6 ```

2. Prevents diffusion into the precursor layer:
   ``` p_2 = \\left(\\frac{h^*}{h + h^*}\\right)^6 ```

The combined pressure for `\\rho_B` is `p = p_1 + p_2`.

# Returns

- `nothing`: Updates `state.rho_B_pressure` with the calculated artificial pressure.
"""
function rho_B_pressure!(state::Active_1D, sys::SysConstActive_1D)
    state.rho_B_pressure .=  (
				  .- Swalbe.power_6.(sys.hcrit ./(state.rho_B .+ sys.b_B))
				  .+ Swalbe.power_6.(sys.hcrit ./(state.height .+ sys.hcrit2))
			)
end

