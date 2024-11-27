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
    @. slipx .= (6μ * height * velx) / (2 * height^2 + 6δ * height + 3δ^2)
    @. slipy .= (6μ * height * vely) / (2 * height^2 + 6δ * height + 3δ^2)
    return nothing
end

slippage!(state::LBM_state_2D, sys::SysConst) = slippage!(
    state.slipx,
    state.slipy,
    state.height,
    state.velx,
    state.vely,
    sys.param.δ,
    sys.param.μ,
)

slippage!(state::Expanded_2D, sys::SysConst) = slippage!(
    state.basestate.slipx,
    state.basestate.slipy,
    state.basestate.height,
    state.basestate.velx,
    state.basestate.vely,
    sys.param.δ,
    sys.param.μ,
)

function slippage!(slip, height, vel, δ, μ)
    @. slip .= (6μ * height * vel) / (2 * height^2 + 6δ * height + 3δ^2)
    return nothing
end

slippage!(state::LBM_state_1D, sys::Consts_1D) = slippage!(state.slip, state.height, state.vel, sys.param.δ, sys.param.μ)

slippage!(state::Expanded_1D, sys::Consts_1D) = slippage!(
    state.basestate.slip,
    state.basestate.height,
    state.basestate.vel,
    sys.param.δ,
    sys.param.μ,
)

slippage!(state::Active_1D, sys::SysConstActive_1D) = slippage!(state.slip, state.height, state.vel, sys.δ, sys.μ)
slippage!(state::LBM_state_1D, sys::SysConstMiscible_1D) = slippage!(state.slip, state.height, state.vel, sys.delta, sys.μ)

# Dirty hack for reducing slip length
function slippage2!(state::LBM_state_2D, sys::SysConst)
    @. state.slipx .=
        (6sys.param.μ * (state.height + sys.param.hcrit) * state.velx) / (
            2 * (state.height + sys.param.hcrit)^2 +
            6sys.param.δ * (state.height + sys.param.hcrit) +
            3sys.param.δ^2
        )
    @. state.slipy .=
        (6sys.param.μ * (state.height + sys.param.hcrit) * state.vely) / (
            2 * (state.height + sys.param.hcrit)^2 +
            6sys.param.δ * (state.height + sys.param.hcrit) +
            3sys.param.δ^2
        )
    return nothing
end



"""
	function slippage!(state::MultiLayer_2D, sys::SysConstMultiLayer)

Calculates the friction force for a two layer thin film according to [Richter et al](https://arxiv.org/abs/2409.16659). There is a maple script generating the used code under `SI/multilayer_symbolics.mv`

Right now no implementation of 3 layers or higher in 3D exists, in 2D we have a three layer implementation. 
"""
function slippage!(state::MultiLayer_2D, sys::SysConstMultiLayer)
    # if sys.layers==2
    #save some combinations we will need multiple times
    state.hi[:,:,1] .= state.height[:,:,1] .+ state.height[:,:,2]
    z_1z_2, z_iz_i = viewneighborsMultiLayer(state.dgrad)
    z_1z_2[:,:,1] .= state.height[:,:,1] .* state.hi[:,:,1] 
    z_iz_i[:,:,1] .= state.height[:,:,1] .* state.height[:,:,1]
    z_iz_i[:,:,2] .= state.hi[:,:,1] .* state.hi[:,:,1]


        state.slipx[:,:, :] .=   (
         1 ./(
               9*sys.delta[1]*sys.delta[1]*sys.delta[2]*sys.delta[2]*sys.mu[1]
               .+ state.height[:,:,1] .* (
                    18*sys.delta[1]*sys.delta[2]*sys.mu[1]*(-sys.delta[1]+sys.delta[2])
               )
               .+ state.hi[:,:,1] .* ( 18*sys.delta[1]*sys.delta[1]*sys.delta[2]*sys.mu[1])
               .+ z_iz_i[:,:,1] .* (
                    6*sys.delta[1]*sys.delta[1]*sys.mu[1]-9*sys.delta[1]*sys.delta[1]*sys.mu[2]-36*sys.delta[1]*sys.delta[2]*sys.mu[1]+6*sys.delta[2]*sys.delta[2]*sys.mu[1]
               )
               .+ z_1z_2[:,:,1] .* (
                    -12*sys.delta[1]*sys.delta[1]*sys.mu[1]+9*sys.delta[1]*sys.delta[1]*sys.mu[2]+36*sys.delta[1]*sys.delta[2]*sys.mu[1]
               )
               .+ z_iz_i[:,:,2] .* (
                    6*sys.delta[1]*sys.delta[1]*sys.mu[1]
               )
               .+ z_iz_i[:,:,1] .* state.height[:,:,1] .* (
                12*(sys.delta[1]*sys.mu[1]-sys.delta[1]*sys.mu[2]-sys.delta[2]*sys.mu[1])
               )
               .+ z_iz_i[:,:,1] .* state.hi[:,:,1] .* (
                    -24*sys.delta[1]*sys.mu[1]+12*sys.delta[1]*sys.mu[2]+12*sys.delta[2]*sys.mu[1]
               )
               .+ z_1z_2[:,:,1] .* state.hi[:,:,1] .* (
                    12*sys.delta[1]*sys.mu[1]
               )
               .+ z_iz_i[:,:,1] .* z_iz_i[:,:,1] .* (
                    4*sys.mu[1]-3*sys.mu[2]
               )
               .+ z_iz_i[:,:,1] .* z_1z_2[:,:,1] .* (
                -8*sys.mu[1]+3*sys.mu[2]
               )
               .+ 4*sys.mu[1] .* z_iz_i[:,:,1] .* z_iz_i[:,:,2]
           )
        )

        state.slipy[:,:, :].= state.slipx[:,:,:]
 
       
        state.slipx[:,:,1] .*=  6*sys.mu[1] .* state.height[:,:,1].* (
              state.velx[:,:,1] .* (
                -3*sys.delta[2]*sys.delta[2]*sys.mu[1]
                .- state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]+6*sys.delta[2]*sys.mu[1]
                )
                .+ z_iz_i[:,:,1] .* (
                    -2*sys.mu[1]+6*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    4*sys.mu[1]-6*sys.mu[2]
                )
                .+ z_iz_i[:,:,2] .* (
                    -2*sys.mu[1]
                )
              )
              .+ state.velx[:,:,2] .* (
                state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]
                )
                .+ z_iz_i[:,:,1] .* (
                    -3*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    3*sys.mu[2]
                )
              )
           ) 
     state.slipx[:,:,2] .*=  6*sys.mu[2]*sys.mu[1] .* state.height[:,:,2] .* (
        state.velx[:,:,1]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 3 .* z_iz_i[:,:,1]  
            ) 
        .- state.velx[:,:,2]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 2 .* z_iz_i[:,:,1]
            ) 
     )


        state.slipy[:,:,1] .*=  6*sys.mu[1] .* state.height[:,:,1].* (
              state.vely[:,:,1] .* (
                -3*sys.delta[2]*sys.delta[2]*sys.mu[1]
                .- state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]+6*sys.delta[2]*sys.mu[1]
                )
                .+ z_iz_i[:,:,1] .* (
                    -2*sys.mu[1]+6*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    4*sys.mu[1]-6*sys.mu[2]
                )
                .+ z_iz_i[:,:,2] .* (
                    -2*sys.mu[1]
                )
              )
              .+ state.vely[:,:,2] .* (
                state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]
                )
                .+ z_iz_i[:,:,1] .* (
                    -3*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    3*sys.mu[2]
                )
              )
           ) 
     state.slipy[:,:,2] .*=  6*sys.mu[2]*sys.mu[1] .* state.height[:,:,2] .* (
        state.vely[:,:,1]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 3 .* z_iz_i[:,:,1]  
            ) 
        .- state.vely[:,:,2]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 2 .* z_iz_i[:,:,1]
            ) 
     )
    # end
    return nothing
end


"""
	function slippage!(state::MultiLayer_2D, sys::SysConstMultiLayer)

Calculates the friction force for a two layer thin film according to [Richter et al](https://arxiv.org/abs/2409.16659). There is a maple script generating the used code under `SI/multilayer_symbolics.mv`

Right now no implementation of 3 layers or higher in 3D exists, in 2D we have a three layer implementation. 
"""
function slippage!(state::MultiLayer_2D, sys::SysConstMultiLayer)
    # if sys.layers==2
    #save some combinations we will need multiple times
    state.hi[:,:,1] .= state.height[:,:,1] .+ state.height[:,:,2]
    z_1z_2, z_iz_i = viewneighborsMultiLayer(state.dgrad)
    z_1z_2[:,:,1] .= state.height[:,:,1] .* state.hi[:,:,1] 
    z_iz_i[:,:,1] .= state.height[:,:,1] .* state.height[:,:,1]
    z_iz_i[:,:,2] .= state.hi[:,:,1] .* state.hi[:,:,1]


        state.slipx[:,:, :] .=   (
         1 ./(
               9*sys.delta[1]*sys.delta[1]*sys.delta[2]*sys.delta[2]*sys.mu[1]
               .+ state.height[:,:,1] .* (
                    18*sys.delta[1]*sys.delta[2]*sys.mu[1]*(-sys.delta[1]+sys.delta[2])
               )
               .+ state.hi[:,:,1] .* ( 18*sys.delta[1]*sys.delta[1]*sys.delta[2]*sys.mu[1])
               .+ z_iz_i[:,:,1] .* (
                    6*sys.delta[1]*sys.delta[1]*sys.mu[1]-9*sys.delta[1]*sys.delta[1]*sys.mu[2]-36*sys.delta[1]*sys.delta[2]*sys.mu[1]+6*sys.delta[2]*sys.delta[2]*sys.mu[1]
               )
               .+ z_1z_2[:,:,1] .* (
                    -12*sys.delta[1]*sys.delta[1]*sys.mu[1]+9*sys.delta[1]*sys.delta[1]*sys.mu[2]+36*sys.delta[1]*sys.delta[2]*sys.mu[1]
               )
               .+ z_iz_i[:,:,2] .* (
                    6*sys.delta[1]*sys.delta[1]*sys.mu[1]
               )
               .+ z_iz_i[:,:,1] .* state.height[:,:,1] .* (
                12*(sys.delta[1]*sys.mu[1]-sys.delta[1]*sys.mu[2]-sys.delta[2]*sys.mu[1])
               )
               .+ z_iz_i[:,:,1] .* state.hi[:,:,1] .* (
                    -24*sys.delta[1]*sys.mu[1]+12*sys.delta[1]*sys.mu[2]+12*sys.delta[2]*sys.mu[1]
               )
               .+ z_1z_2[:,:,1] .* state.hi[:,:,1] .* (
                    12*sys.delta[1]*sys.mu[1]
               )
               .+ z_iz_i[:,:,1] .* z_iz_i[:,:,1] .* (
                    4*sys.mu[1]-3*sys.mu[2]
               )
               .+ z_iz_i[:,:,1] .* z_1z_2[:,:,1] .* (
                -8*sys.mu[1]+3*sys.mu[2]
               )
               .+ 4*sys.mu[1] .* z_iz_i[:,:,1] .* z_iz_i[:,:,2]
           )
        )

        state.slipy[:,:, :].= state.slipx[:,:,:]
 
       
        state.slipx[:,:,1] .*=  6*sys.mu[1] .* state.height[:,:,1].* (
              state.velx[:,:,1] .* (
                -3*sys.delta[2]*sys.delta[2]*sys.mu[1]
                .- state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]+6*sys.delta[2]*sys.mu[1]
                )
                .+ z_iz_i[:,:,1] .* (
                    -2*sys.mu[1]+6*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    4*sys.mu[1]-6*sys.mu[2]
                )
                .+ z_iz_i[:,:,2] .* (
                    -2*sys.mu[1]
                )
              )
              .+ state.velx[:,:,2] .* (
                state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]
                )
                .+ z_iz_i[:,:,1] .* (
                    -3*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    3*sys.mu[2]
                )
              )
           ) 
     state.slipx[:,:,2] .*=  6*sys.mu[2]*sys.mu[1] .* state.height[:,:,2] .* (
        state.velx[:,:,1]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 3 .* z_iz_i[:,:,1]  
            ) 
        .- state.velx[:,:,2]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 2 .* z_iz_i[:,:,1]
            ) 
     )


        state.slipy[:,:,1] .*=  6*sys.mu[1] .* state.height[:,:,1].* (
              state.vely[:,:,1] .* (
                -3*sys.delta[2]*sys.delta[2]*sys.mu[1]
                .- state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]+6*sys.delta[2]*sys.mu[1]
                )
                .+ z_iz_i[:,:,1] .* (
                    -2*sys.mu[1]+6*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    4*sys.mu[1]-6*sys.mu[2]
                )
                .+ z_iz_i[:,:,2] .* (
                    -2*sys.mu[1]
                )
              )
              .+ state.vely[:,:,2] .* (
                state.height[:,:,2] .* (
                    6*sys.delta[1]*sys.mu[2]
                )
                .+ z_iz_i[:,:,1] .* (
                    -3*sys.mu[2]
                )
                .+ z_1z_2[:,:,1] .* (
                    3*sys.mu[2]
                )
              )
           ) 
     state.slipy[:,:,2] .*=  6*sys.mu[2]*sys.mu[1] .* state.height[:,:,2] .* (
        state.vely[:,:,1]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 3 .* z_iz_i[:,:,1]  
            ) 
        .- state.vely[:,:,2]  .* (
                3*sys.delta[1]*sys.delta[1] .+ 6*sys.delta[1] .* state.height[:,:,1] .+ 2 .* z_iz_i[:,:,1]
            ) 
     )
    # end
    return nothing
end



"""
	function slippage_no_slip!(state::MultiLayer_2D, sys::SysConstMultiLayer)

Calculates the friction force for a two layer thin film according to [Richter et al](https://arxiv.org/abs/2409.16659) in simplified no-slip case delta=0. There is a maple script generating the used code under `SI/multilayer_symbolics.mv`

Right now no implementation of 3 layers or higher in 3D exists, in 2D we have a three layer implementation. 
"""

function slippage_no_slip!(state::MultiLayer_2D, sys::SysConstMultiLayer)
    # if sys.layers==2
    state.hi[:,:,1] .= state.height[:,:,1] .+ state.height[:,:,2]
        state.slipx[:,:,1] .=  2*sys.mu[1] .* state.height[:,:,1] .* (
            (
                -3 .* (state.velx[:,:,1] .*(state.height[:,:,1] .* (2*sys.mu[1] -6*sys.mu[2]) .- 2*sys.mu[1] .* state.hi[:,:,1] ) .+ 3* sys.mu[2] .* state.velx[:,:,2] .* state.height[:,:,1])
            )./(
                state.height[:,:,1] .* state.height[:,:,1] .* ( state.height[:,:,1] .* (4*sys.mu[1]-3*sys.mu[2]) .- 4*sys.mu[1] .* state.hi[:,:,1])
            ) 
        )
     state.slip[:,:,2] .= 2*sys.mu[2] .* state.height[:,:,2] .* (
        (
            3*sys.mu[1] .* (3 .*state.velx[:,:,1] .- 2 .* state.velx[:,:,2]) 
        )./(
           state.height[:,:,1] .* state.height[:,:,1] .* (4*sys.mu[1]-3*sys.mu[2]) .+ state.height[:,:,1] .* state.hi[:,:,1] .* ( -8*sys.mu[1] + 3*sys.mu[2]) .+ 4*sys.mu[1] .* state.hi[:,:,1] .* state.hi[:,:,1] 
        )
     )
     state.slipy[:,:,1] .=  2*sys.mu[1] .* state.height[:,:,1] .* (
            (
                -3 .* (state.vely[:,:,1] .*(state.height[:,:,1] .* (2*sys.mu[1] -6*sys.mu[2]) .- 2*sys.mu[1] .* state.hi[:,:,1] ) .+ 3* sys.mu[2] .* state.vely[:,:,2] .* state.height[:,:,1])
            )./(
                state.height[:,:,1] .* state.height[:,:,1] .* ( state.height[:,:,1] .* (4*sys.mu[1]-3*sys.mu[2]) .- 4*sys.mu[1] .* state.hi[:,:,1])
            ) 
        )
     state.slip[:,:,2] .= 2*sys.mu[2] .* state.height[:,:,2] .* (
        (
            3*sys.mu[1] .* (3 .*state.vely[:,:,1] .- 2 .* state.vely[:,:,2]) 
        )./(
           state.height[:,:,1] .* state.height[:,:,1] .* (4*sys.mu[1]-3*sys.mu[2]) .+ state.height[:,:,1] .* state.hi[:,:,1] .* ( -8*sys.mu[1] + 3*sys.mu[2]) .+ 4*sys.mu[1] .* state.hi[:,:,1] .* state.hi[:,:,1] 
        )
     )
    # end
    return nothing
end

"""
    function slippage!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)

Calculates the friction force for a multilayer system with slip.

    # Variables
`state :: StateMultiLayer_1D` The fields
`sys :: SysConstMultiLayer_1D` The system Variables

    # Mathematics

    ## 2 layers
    The friction force reads ``F_i=2h_i\\mu_i a_i``

    kwhere ``a_i`` is the solution of the linear system
    ```math
      \\begin{pmatrix}
    \\beta_1^2 & - \\beta_1 & 1 & 0 & 0 & 0 \\\\
    \\frac{z_1^2}{3} & \\frac{z_1}{2} & 1 & 0 & 0 & 0 \\\\
    2 \\mu_1 z_1 & \\mu_1 & 0 & -2 \\mu_2 z_1 & - \\mu_2 & 0 \\\\
    - z_1^2 & - z_1 & -1 & (z_1-\\beta_2)^2  & z_1-\\beta_2 & 1 \\\\
    0& 0& 0&\\frac{z_2^3-z_1^3}{3(z_2-z_1)} & \\frac{z_2^2-z_1^2}{2(z_2-z_1)} & 1\\\\
    0 & 0 & 0 & 2z_2 & 1 & 0  
    \\end{pmatrix}\\begin{pmatrix}
    a_1 \\\\ b_1 \\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2
    \\end{pmatrix}=\\begin{pmatrix}
    0 \\\\ U_1 \\\\ 0 \\\\ 0 \\\\  U_2\\\\ 0
    \\end{pmatrix}. 
    ```

    ## 3 layers 
    Actually solving the linear system gives ridiculously long terms. We do the next best thing, solving 
    ```math
    \\begin{pmatrix}
        0 e 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\\\
        \\frac{z_1^3}{3} & \\frac{z_1^2}{2} & z_1 & 0 & 0 & 0& 0 & 0 & 0\\\\
        2z_1 & 1 & 0 & -2\\frac{\\mu_2}{\\mu_1}z_1 & - \\frac{\\mu_2}{\\mu_1 } & 0 & 0 & 0 & 0\\\\
        z_1^2 & z_1 & 1 & - z_1^2 &- z_1 & -1 & 0 & 0 & 0\\\\
        0 & 0 & 0 & \\frac{z_2^3-z_1^3}{3 }& \\frac{z_2^2-z_1^2}{2} & z_2-z_1 & 0 & 0 & 0\\\\
        0 & 0 & 0 & 2z_2 & 1 & 0 & -2\\frac{\\mu_3}{\\mu_2} z_2 & -\\frac{\\mu_3}{\\mu_2} & 0 \\\\
        0 & 0 & 0 & z_2^2 & z_2 & 1 & -z_2^2 & - z_2 & -1 \\\\
        0 & 0 & 0 & 0 & 0 & 0 & \\frac{z_3^3-z_2^3}{3} & \\frac{z_3^2 - z_2^3}{2} & z_3-z_2\\\\
        0 & 0 & 0 & 0 & 0 & 0 & 2z_3 & 1 & 0 
    \\end{pmatrix}
    \\begin{pmatrix}
        a_1 \\\\ b_1 \\\\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2 \\\\ a_3 \\\\ b_3 \\\\ c_3 
    \\end{pmatrix}
    = \\begin{pmatrix}
        0 \\\\ z_1 U_1 \\\\ 0 \\\\ 0 \\\\ (z_2-z_1) U_2 \\\\ 0 \\\\ 0\\\\ (z_3-z_2) U_3 \\\\ 0
    \\end{pmatrix}
    ```
    and then we just add a factor of ``\\mu_i\\mu_jd^n`` under the fraction such that the dimensions are correct. That is not the actual thing but just a very pragmatic solution. Whenever you are really interested in the dynamics You better resort to `slippage_no_slip!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)`

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.

# Number of layers

With slip we have only implemented two layers, there is a three layer implementation of the simplified no-slip case [slippage_no_slip!](@ref)
"""
function slippage!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)
    if sys.layers==2
        state.hi[:,1] .= state.height[:,1] .+ state.height[:,2]
        state.slip[:,1] .= 2*sys.mu[1] .* state.height[:,1] .* friction_force_2_1_slip.(sys.mu[1],sys.mu[2], state.height[:,1], state.hi[:,1], state.vel[:,1], state.vel[:,2],sys.delta[1],sys.delta[2])
        state.slip[:,2] .= 2*sys.mu[2] .* state.height[:,2] .* friction_force_2_2_slip.(sys.mu[1],sys.mu[2], state.height[:,1], state.hi[:,1], state.vel[:,1], state.vel[:,2],sys.delta[1],sys.delta[2])
    elseif sys.layers==3
        state.hi[:,1] .= state.height[:,1] .+ state.height[:,2]
        state.hi[:,2] .= state.hi[:,1] .+ state.height[:,3]
        state.slip[:,1] .= 2*sys.mu[1] .* state.height[:,1] .* friction_force_3_1_slip.(sys.mu[1], sys.mu[2], sys.mu[3], state.height[:,1], state.hi[:,1], state.hi[:,2], state.vel[:,1], state.vel[:,2], state.vel[:,3],sys.delta[1])
        state.slip[:,2] .= 2*sys.mu[2] .* state.height[:,2] .* friction_force_3_2_slip.(sys.mu[1], sys.mu[2], sys.mu[3], state.height[:,1], state.hi[:,1], state.hi[:,2], state.vel[:,1], state.vel[:,2], state.vel[:,3],sys.delta[1])
        state.slip[:,3] .= 2*sys.mu[3] .* state.height[:,3] .* friction_force_3_3_slip.(sys.mu[1], sys.mu[2], sys.mu[3], state.height[:,1], state.hi[:,1], state.hi[:,2], state.vel[:,1], state.vel[:,2], state.vel[:,3],sys.delta[1])
    end
    return nothing
end




"""
    function slippage_no_slip!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)

Calculates the friction force for a multilayer system with slip.

    # Variables
`state :: StateMultiLayer_1D` The fields
`sys :: SysConstMultiLayer_1D` The system Variables

    # Mathematics

    ## 2 layers
    The friction force reads ``F_i=2h_i\\mu_i a_i``

    where ``a_i`` is the solution of the linear system
    ```math
      \\begin{pmatrix}
    0 & 0 & 1 & 0 & 0 & 0 \\\\
    \\frac{z_1^2}{3} & \\frac{z_1}{2} & 1 & 0 & 0 & 0 \\\\
    2 \\mu_1 z_1 & \\mu_1 & 0 & -2 \\mu_2 z_1 & - \\mu_2 & 0 \\\\
    - z_1^2 & - z_1 & -1 & z_1^2  & z_1& 1 \\\\
    0& 0& 0&\\frac{z_2^3-z_1^3}{3(z_2-z_1)} & \\frac{z_2^2-z_1^2}{2(z_2-z_1)} & 1\\\\
    0 & 0 & 0 & 2z_2 & 1 & 0  
    \\end{pmatrix}\\begin{pmatrix}
    a_1 \\\\ b_1 \\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2
    \\end{pmatrix}=\\begin{pmatrix}
    0 \\\\ U_1 \\\\ 0 \\\\ 0 \\\\  U_2\\\\ 0
    \\end{pmatrix}. 
    ```

    ## 3 layers 
    Actually solving the linear system gives ridiculously long terms. We do the next best thing, solving 
    ```math
    \\begin{pmatrix}
        0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\\\
        \\frac{z_1^3}{3} & \\frac{z_1^2}{2} & z_1 & 0 & 0 & 0& 0 & 0 & 0\\\\
        2z_1 & 1 & 0 & -2\\frac{\\mu_2}{\\mu_1}z_1 & - \\frac{\\mu_2}{\\mu_1 } & 0 & 0 & 0 & 0\\\\
        z_1^2 & z_1 & 1 & - z_1^2 &- z_1 & -1 & 0 & 0 & 0\\\\
        0 & 0 & 0 & \\frac{z_2^3-z_1^3}{3 }& \\frac{z_2^2-z_1^2}{2} & z_2-z_1 & 0 & 0 & 0\\\\
        0 & 0 & 0 & 2z_2 & 1 & 0 & -2\\frac{\\mu_3}{\\mu_2} z_2 & -\\frac{\\mu_3}{\\mu_2} & 0 \\\\
        0 & 0 & 0 & z_2^2 & z_2 & 1 & -z_2^2 & - z_2 & -1 \\\\
        0 & 0 & 0 & 0 & 0 & 0 & \\frac{z_3^3-z_2^3}{3} & \\frac{z_3^2 - z_2^3}{2} & z_3-z_2\\\\
        0 & 0 & 0 & 0 & 0 & 0 & 2z_3 & 1 & 0 
    \\end{pmatrix}
    \\begin{pmatrix}
        a_1 \\\\ b_1 \\\\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2 \\\\ a_3 \\\\ b_3 \\\\ c_3 
    \\end{pmatrix}
    = \\begin{pmatrix}
        0 \\\\ z_1 U_1 \\\\ 0 \\\\ 0 \\\\ (z_2-z_1) U_2 \\\\ 0 \\\\ 0\\\\ (z_3-z_2) U_3 \\\\ 0
    \\end{pmatrix}
    ```
    and then we just add a factor of ``\\mu_i\\mu_jd^n`` under the fraction such that the dimensions are correct. That is not the actual thing but just a very pragmatic solution. Whenever you are really interested in the dynamics You better resort to `slippage_no_slip!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)`

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.

# Number of layers

Here we have implemented also the three layer case

"""
function slippage_no_slip!(state::StateMultiLayer_1D, sys::SysConstMultiLayer_1D)
    if sys.layers==2
       state.hi[:,1] .= state.height[:,1] .+ state.height[:,2]
        state.slip[:,1] .= 2*sys.mu[1] .* state.height[:,1] .* friction_force_2_1_no_slip.(sys.mu[1],sys.mu[2], state.height[:,1], state.hi[:,1], state.vel[:,1], state.vel[:,2])
        state.slip[:,2] .= 2*sys.mu[2] .* state.height[:,2] .* friction_force_2_2_no_slip.(sys.mu[1],sys.mu[2], state.height[:,1], state.hi[:,1], state.vel[:,1], state.vel[:,2])
    elseif sys.layers==3
        state.hi[:,1] .= state.height[:,1] .+ state.height[:,2]
        state.hi[:,2] .= state.hi[:,1] .+ state.height[:,3]
        state.slip[:,1] .= 2*sys.mu[1] .* state.height[:,1] .* friction_force_3_1_no_slip.(sys.mu[1], sys.mu[2], sys.mu[3], state.height[:,1], state.hi[:,1], state.hi[:,2], state.vel[:,1], state.vel[:,2], state.vel[:,3])
        state.slip[:,2] .= 2*sys.mu[2] .* state.height[:,2] .* friction_force_3_2_no_slip.(sys.mu[1], sys.mu[2], sys.mu[3], state.height[:,1], state.hi[:,1], state.hi[:,2], state.vel[:,1], state.vel[:,2], state.vel[:,3])
        state.slip[:,3] .= 2*sys.mu[3] .* state.height[:,3] .* friction_force_3_3_no_slip.(sys.mu[1], sys.mu[2], sys.mu[3], state.height[:,1], state.hi[:,1], state.hi[:,2], state.vel[:,1], state.vel[:,2], state.vel[:,3])
    end
    return nothing
end

"""
    function friction_force_2_1_no_slip(mu1,mu2,z1,z2,U1,U2)

Returns the coefficient a_1 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

    # Arguments

-`mu1 :: Float64` first viscosity  
-`mu2 :: Float64` second viscosity  
-`z1 :: Float64` first interface position 
-`z2 :: Float64` second interface position 
-`U1 :: Float64` first flux
-`U2 :: Float64` second flux

    # Mathematics

This is the frirst entry of the solution of 

```math
\\begin{pmatrix}
    0 & 0 & 1 & 0 & 0 & 0 \\\\
    \\frac{z_1^2}{3} & \\frac{z_1}{2} & 1 & 0 & 0 & 0 \\\\
    2 \\mu_1 z_1 & \\mu_1 & 0 & -2 \\mu_2 z_1 & - \\mu_2 & 0 \\\\
    - z_1^2 & - z_1 & -1 & z_1^2  & z_1 & 1 \\\\
    0& 0& 0&\\frac{z_2^3-z_1^3}{3(z_2-z_1)} & \\frac{z_2^2-z_1^2}{2(z_2-z_1)} & 1\\\\
    0 & 0 & 0 & 2z_2 & 1 & 0  
    \\end{pmatrix}\\begin{pmatrix}
    a_1 \\\\ b_1 \\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2
    \\end{pmatrix}=\\begin{pmatrix}
    0 \\\\ U_1 \\\\ 0 \\\\ 0 \\\\  U_2\\\\ 0
    \\end{pmatrix}. 
```

    # Comments

The linear system has been solved by maple and copy pasted here. 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_2_1_no_slip(mu1,mu2,z1,z2,U1,U2)
    return   -6 * (mu1 * z1 - mu1 * z2 - 3 * mu2 * z1) / (4 * mu1 * z1 - 4 * mu1 * z2 - 3 * mu2 * z1) / z1 ^ 2 * U1 - 9 / z1 * mu2 / (4 * mu1 * z1 - 4 * mu1 * z2 - 3 * mu2 * z1) * U2
end
"""
    function friction_force_2_2_no_slip(mu1,mu2,z1,z2,U1,U2)

Returns the coefficient a_1 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

  see `friction_force_2_1_no_slip(mu1,mu2,z1,z2,U1,U2)`
```

    # Comments

The linear system has been solved by maple and copy pasted here. 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_2_2_no_slip(mu1,mu2,z1,z2,U1,U2)
    return   9 * mu1 / (4 * z1 ^ 2 * mu1 - 8 * z2 * z1 * mu1 + 4 * mu1 * z2 ^ 2 - 3 * z1 ^ 2 * mu2 + 3 * z2 * z1 * mu2) * U1 - 6 * mu1 / (4 * z1 ^ 2 * mu1 - 8 * z2 * z1 * mu1 + 4 * mu1 * z2 ^ 2 - 3 * z1 ^ 2 * mu2 + 3 * z2 * z1 * mu2) * U2
end


"""
    function friction_force_2_1_slip(mu1,mu2,z1,z2,U1,U2)

Returns the coefficient a_1 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

    # Arguments

-`mu1 :: Float64` first viscosity  
-`mu2 :: Float64` second viscosity  
-`z1 :: Float64` first interface position 
-`z2 :: Float64` second interface position 
-`U1 :: Float64` first flux
-`U2 :: Float64` second flux
-`d1 :: Float64` first slip length
-`d2 :: Float64` second slip length

    # Mathematics

This is the frirst entry of the solution of 

```math
\\beta_1^2 & - \\beta_1 & 1 & 0 & 0 & 0 \\\\
\\frac{z_1^2}{3} & \\frac{z_1}{2} & 1 & 0 & 0 & 0 \\\\
2 \\mu_1 z_1 & \\mu_1 & 0 & -2 \\mu_2 z_1 & - \\mu_2 & 0 \\\\
- z_1^2 & - z_1 & -1 & (z_1-\\beta_2)^2  & z_1-\\beta_2 & 1 \\\\
0& 0& 0&\\frac{z_2^3-z_1^3}{3(z_2-z_1)} & \\frac{z_2^2-z_1^2}{2(z_2-z_1)} & 1\\\\
0 & 0 & 0 & 2z_2 & 1 & 0  
\\end{pmatrix}\\begin{pmatrix}
a_1 \\\\ b_1 \\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2
\\end{pmatrix}=\\begin{pmatrix}
0 \\\\ U_1 \\\\ 0 \\\\ 0 \\\\  U_2\\\\ 0
\\end{pmatrix}.
```

    # Comments

The linear system has been solved by maple and copy pasted here. 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_2_1_slip(mu1,mu2,z1,z2,U1,U2,d1,d2)
    return  (
        3 * (6 * d1 * mu2 * z1 - 6 * d1 * mu2 * z2 - 3 * d2 ^ 2 * mu1 + 6 * d2 * mu1 * z1 - 6 * d2 * mu1 * z2 - 2 * z1 ^ 2 * mu1 + 4 * z2 * z1 * mu1 - 2 * mu1 * z2 ^ 2 + 6 * z1 ^ 2 * mu2 - 6 * z2 * z1 * mu2) / (9 * d1 ^ 2 * d2 ^ 2 * mu1 - 18 * d1 ^ 2 * d2 * mu1 * z1 + 18 * d1 ^ 2 * d2 * mu1 * z2 + 6 * d1 ^ 2 * mu1 * z1 ^ 2 - 12 * d1 ^ 2 * mu1 * z1 * z2 + 6 * d1 ^ 2 * mu1 * z2 ^ 2 - 9 * d1 ^ 2 * mu2 * z1 ^ 2 + 9 * d1 ^ 2 * mu2 * z1 * z2 + 18 * d1 * d2 ^ 2 * mu1 * z1 - 36 * d1 * d2 * mu1 * z1 ^ 2 + 36 * d1 * d2 * mu1 * z1 * z2 + 12 * d1 * mu1 * z1 ^ 3 - 24 * d1 * mu1 * z1 ^ 2 * z2 + 12 * d1 * mu1 * z1 * z2 ^ 2 - 12 * d1 * mu2 * z1 ^ 3 + 12 * d1 * mu2 * z1 ^ 2 * z2 + 6 * d2 ^ 2 * mu1 * z1 ^ 2 - 12 * d2 * mu1 * z1 ^ 3 + 12 * d2 * mu1 * z1 ^ 2 * z2 + 4 * mu1 * z1 ^ 4 - 8 * mu1 * z1 ^ 3 * z2 + 4 * mu1 * z1 ^ 2 * z2 ^ 2 - 3 * mu2 * z1 ^ 4 + 3 * mu2 * z1 ^ 3 * z2) * U1
         - 9 * mu2 * (z1 - z2) * (2 * d1 + z1) / (9 * d1 ^ 2 * d2 ^ 2 * mu1 - 18 * d1 ^ 2 * d2 * mu1 * z1 + 18 * d1 ^ 2 * d2 * mu1 * z2 + 6 * d1 ^ 2 * mu1 * z1 ^ 2 - 12 * d1 ^ 2 * mu1 * z1 * z2 + 6 * d1 ^ 2 * mu1 * z2 ^ 2 - 9 * d1 ^ 2 * mu2 * z1 ^ 2 + 9 * d1 ^ 2 * mu2 * z1 * z2 + 18 * d1 * d2 ^ 2 * mu1 * z1 - 36 * d1 * d2 * mu1 * z1 ^ 2 + 36 * d1 * d2 * mu1 * z1 * z2 + 12 * d1 * mu1 * z1 ^ 3 - 24 * d1 * mu1 * z1 ^ 2 * z2 + 12 * d1 * mu1 * z1 * z2 ^ 2 - 12 * d1 * mu2 * z1 ^ 3 + 12 * d1 * mu2 * z1 ^ 2 * z2 + 6 * d2 ^ 2 * mu1 * z1 ^ 2 - 12 * d2 * mu1 * z1 ^ 3 + 12 * d2 * mu1 * z1 ^ 2 * z2 + 4 * mu1 * z1 ^ 4 - 8 * mu1 * z1 ^ 3 * z2 + 4 * mu1 * z1 ^ 2 * z2 ^ 2 - 3 * mu2 * z1 ^ 4 + 3 * mu2 * z1 ^ 3 * z2) * U2
    )
end
"""
    function friction_force_2_2_slip(mu1,mu2,z1,z2,U1,U2)

Returns the coefficient a_1 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

  see `friction_force_2_1_slip(mu1,mu2,z1,z2,U1,U2)`
```

    # Comments

The linear system has been solved by maple and copy pasted here. 


# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_2_2_slip(mu1,mu2,z1,z2,U1,U2,d1,d2)
    return  (
        9 * mu1 * (d1 ^ 2 + 2 * d1 * z1 + z1 ^ 2) / (9 * d1 ^ 2 * d2 ^ 2 * mu1 - 18 * d1 ^ 2 * d2 * mu1 * z1 + 18 * d1 ^ 2 * d2 * mu1 * z2 + 6 * d1 ^ 2 * mu1 * z1 ^ 2 - 12 * d1 ^ 2 * mu1 * z1 * z2 + 6 * d1 ^ 2 * mu1 * z2 ^ 2 - 9 * d1 ^ 2 * mu2 * z1 ^ 2 + 9 * d1 ^ 2 * mu2 * z1 * z2 + 18 * d1 * d2 ^ 2 * mu1 * z1 - 36 * d1 * d2 * mu1 * z1 ^ 2 + 36 * d1 * d2 * mu1 * z1 * z2 + 12 * d1 * mu1 * z1 ^ 3 - 24 * d1 * mu1 * z1 ^ 2 * z2 + 12 * d1 * mu1 * z1 * z2 ^ 2 - 12 * d1 * mu2 * z1 ^ 3 + 12 * d1 * mu2 * z1 ^ 2 * z2 + 6 * d2 ^ 2 * mu1 * z1 ^ 2 - 12 * d2 * mu1 * z1 ^ 3 + 12 * d2 * mu1 * z1 ^ 2 * z2 + 4 * mu1 * z1 ^ 4 - 8 * mu1 * z1 ^ 3 * z2 + 4 * mu1 * z1 ^ 2 * z2 ^ 2 - 3 * mu2 * z1 ^ 4 + 3 * mu2 * z1 ^ 3 * z2) * U1
         - 3 * (3 * d1 ^ 2 + 6 * d1 * z1 + 2 * z1 ^ 2) * mu1 / (9 * d1 ^ 2 * d2 ^ 2 * mu1 - 18 * d1 ^ 2 * d2 * mu1 * z1 + 18 * d1 ^ 2 * d2 * mu1 * z2 + 6 * d1 ^ 2 * mu1 * z1 ^ 2 - 12 * d1 ^ 2 * mu1 * z1 * z2 + 6 * d1 ^ 2 * mu1 * z2 ^ 2 - 9 * d1 ^ 2 * mu2 * z1 ^ 2 + 9 * d1 ^ 2 * mu2 * z1 * z2 + 18 * d1 * d2 ^ 2 * mu1 * z1 - 36 * d1 * d2 * mu1 * z1 ^ 2 + 36 * d1 * d2 * mu1 * z1 * z2 + 12 * d1 * mu1 * z1 ^ 3 - 24 * d1 * mu1 * z1 ^ 2 * z2 + 12 * d1 * mu1 * z1 * z2 ^ 2 - 12 * d1 * mu2 * z1 ^ 3 + 12 * d1 * mu2 * z1 ^ 2 * z2 + 6 * d2 ^ 2 * mu1 * z1 ^ 2 - 12 * d2 * mu1 * z1 ^ 3 + 12 * d2 * mu1 * z1 ^ 2 * z2 + 4 * mu1 * z1 ^ 4 - 8 * mu1 * z1 ^ 3 * z2 + 4 * mu1 * z1 ^ 2 * z2 ^ 2 - 3 * mu2 * z1 ^ 4 + 3 * mu2 * z1 ^ 3 * z2) * U2
    )
end

"""
    function friction_force_3_1_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)

Returns the coefficient a_1 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

    # Arguments

-`mu1 :: Float64` first viscosity  
-`mu2 :: Float64` second viscosity  
-`mu3 :: Float64` third viscosity  
-`z1 :: Float64` first interface position 
-`z2 :: Float64` second interface position 
-`z3 :: Float64` third interface position 
-`U1 :: Float64` first flux
-`U2 :: Float64` second flux
-`U3 :: Float64` third flux

    # Mathematics

This is the frirst entry of the solution of 

```math
\\begin{pmatrix}
        0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\\\
        \\frac{z_1^3}{3} & \\frac{z_1^2}{2} & z_1 & 0 & 0 & 0& 0 & 0 & 0\\\\
        2z_1 & 1 & 0 & -2\\frac{\\mu_2}{\\mu_1}z_1 & - \\frac{\\mu_2}{\\mu_1 } & 0 & 0 & 0 & 0\\\\
        z_1^2 & z_1 & 1 & - z_1^2 &- z_1 & -1 & 0 & 0 & 0\\\\
        0 & 0 & 0 & \\frac{z_2^3-z_1^3}{3 }& \\frac{z_2^2-z_1^2}{2} & z_2-z_1 & 0 & 0 & 0\\\\
        0 & 0 & 0 & 2z_2 & 1 & 0 & -2\\frac{\\mu_3}{\\mu_2} z_2 & -\\frac{\\mu_3}{\\mu_2} & 0 \\\\
        0 & 0 & 0 & z_2^2 & z_2 & 1 & -z_2^2 & - z_2 & -1 \\\\
        0 & 0 & 0 & 0 & 0 & 0 & \\frac{z_3^3-z_2^3}{3} & \\frac{z_3^2 - z_2^3}{2} & z_3-z_2\\\\
        0 & 0 & 0 & 0 & 0 & 0 & 2z_3 & 1 & 0 
    \\end{pmatrix}
    \\begin{pmatrix}
        a_1 \\\\ b_1 \\\\\\ c_1 \\\\ a_2 \\\\ b_2 \\\\ c_2 \\\\ a_3 \\\\ b_3 \\\\ c_3 
    \\end{pmatrix}
    = \\begin{pmatrix}
        0 \\\\ z_1 U_1 \\\\ 0 \\\\ 0 \\\\ (z_2-z_1) U_2 \\\\ 0 \\\\ 0\\\\ (z_3-z_2) U_3 \\\\ 0
    \\end{pmatrix}
```

    # Comments

The linear system has been solved by maple and copy pasted here. 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_3_1_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)
    return (
        -3/2 * (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 12 * mu2 ^ 2 * z1 * z2 + 12 * mu2 ^ 2 * z1 * z3 - 12 * mu2 * mu3 * z1 ^ 2 + 12 * mu2 * mu3 * z1 * z2) / (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) / z1 ^ 2 * U1
     + 9/2 * mu2 / z1 * (2 * mu2 * z2 - 2 * mu2 * z3 + 3 * mu3 * z1 - 3 * mu3 * z2) / (-z2 + z1) / (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) * (z2 - z1) * U2 
    - 9/2 / z1 * (-z2 + z1) * mu3 * mu2 / (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) * (z3 - z2) * U3
    )
end


"""
    function friction_force_3_2_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)

Returns the coefficient a_2 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

See `friction_force_3_1_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)`

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_3_2_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)
    return (
        9/2 * mu1 * (2 * mu2 * z2 - 2 * mu2 * z3 + 3 * mu3 * z1 - 3 * mu3 * z2) / (-z2 + z1) / (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) * U1
         + 3/2 * (4 * mu1 * mu2 * z2 - 4 * mu1 * mu2 * z3 + 12 * mu1 * mu3 * z1 - 12 * mu1 * mu3 * z2 - 3 * mu2 * mu3 * z1) / (-z2 + z1) ^ 2 / (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) * (z2 - z1) * U2
          - 9/2 * (2 * mu1 * z1 - 2 * mu1 * z2 - mu2 * z1) * mu3 / (-z2 + z1) / (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) * (z3 - z2) * U3
    )
end

"""
    function friction_force_3_2_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)

Returns the coefficient a_3 for the friction force of the lowest (counting from the botton) layer of a three layer system without slip

See `friction_force_3_1_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)` 



# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_3_3_no_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)
   return (
        -9/2 * mu2 * (-z2 + z1) * mu1 / (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) * U1
         - 9/2 * mu2 * (2 * mu1 * z1 - 2 * mu1 * z2 - mu2 * z1) / (-z2 + z1) / (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) * (z2 - z1) * U2
          + 3/2 * (4 * mu1 * z1 - 4 * mu1 * z2 - 3 * mu2 * z1) * mu2 / (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) / (-z3 + z2) * (z3 - z2) * U3
    )
end


"""
    function friction_force_3_1_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3,d)

Returns the coefficient a_1 for the friction force of the lowest (counting from the botton) layer of a three layer system with slip

    # Arguments

-`mu1 :: Float64` first viscosity  
-`mu2 :: Float64` second viscosity  
-`mu3 :: Float64` third viscosity  
-`z1 :: Float64` first interface position 
-`z2 :: Float64` second interface position 
-`z3 :: Float64` third interface position 
-`U1 :: Float64` first flux
-`U2 :: Float64` second flux
-`U3 :: Float64` third flux
-`d :: Float64` slip length

    # Mathematics

Solving the actual equations for a 3 layer film with slippage leads to a ridicusly long term. Instead I just add a term proportional to ``\\mu_i\\mu_j*d^4`` under the fraction 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_3_1_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3,d)
    return (
        -3/2 * (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 12 * mu2 ^ 2 * z1 * z2 + 12 * mu2 ^ 2 * z1 * z3 - 12 * mu2 * mu3 * z1 ^ 2 + 12 * mu2 * mu3 * z1 * z2) / ((4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) * z1 ^ 2 + 1.5 * mu1 * mu1* d^4) * U1
     + 9/2 * mu2 * (2 * mu2 * z2 - 2 * mu2 * z3 + 3 * mu3 * z1 - 3 * mu3 * z2) / (z1 * (-z2 + z1) * (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) + 4.5 * mu1 * mu2 * d^4) * (z2 - z1) * U2 
    - 9/2 * (-z2 + z1) * mu3 * mu2 / ( z1 * (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) + 4.5 * mu1 * mu3* d^4) * (z3 - z2) * U3
    )
end


"""
    function friction_force_3_2_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3,d)

Returns the coefficient a_2 for the friction force of the lowest (counting from the botton) layer of a three layer system with slip

See `friction_force_3_1_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)`

Solving the actual equations for a 3 layer film with slippage leads to a ridicusly long term. Instead I just add a term proportional to ``\\mu_i\\mu_jd^n``` under the fraction where n is chosen such that dimensions are respected 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_3_2_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3,d)
    return (
        9/2 * mu1 * (2 * mu2 * z2 - 2 * mu2 * z3 + 3 * mu3 * z1 - 3 * mu3 * z2)  / ((-z2 + z1) *(4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) + 4.5* mu2 * mu1 *d^3) * U1
         + 3/2 * (4 * mu1 * mu2 * z2 - 4 * mu1 * mu2 * z3 + 12 * mu1 * mu3 * z1 - 12 * mu1 * mu3 * z2 - 3 * mu2 * mu3 * z1) /((-z2 + z1)^2 * (4 * mu1 * mu2 * z1 * z2 - 4 * mu1 * mu2 * z1 * z3 - 4 * mu1 * mu2 * z2 ^ 2 + 4 * mu1 * mu2 * z2 * z3 + 3 * mu1 * mu3 * z1 ^ 2 - 6 * mu1 * mu3 * z1 * z2 + 3 * mu1 * mu3 * z2 ^ 2 - 3 * mu2 ^ 2 * z1 * z2 + 3 * mu2 ^ 2 * z1 * z3 - 3 * mu2 * mu3 * z1 ^ 2 + 3 * mu2 * mu3 * z1 * z2) + 1.5 * mu2 * mu2 * d^4) * (z2 - z1) * U2
          - 9/2 * (2 * mu1 * z1 - 2 * mu1 * z2 - mu2 * z1) * mu3 /((-z2 + z1) * (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) + 4.5 * mu2 * mu3* d^4) * (z3 - z2) * U3
    )
end

"""
    function friction_force_3_2_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3,d)

Returns the coefficient a_3 for the friction force of the lowest (counting from the botton) layer of a three layer system with slip

See `friction_force_3_1_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3)` 

Solving the actual equations for a 3 layer film with slippage leads to a ridicusly long term. Instead I just add a term proportional to ``\\mu_i\\mu_jd^n``` under the fraction where n is chosen such that dimensions are respected 

# Reference

There is a symbolic calculation file creating the here used code at `SI/multilayer_symbolics.mv`.


"""
function friction_force_3_3_slip(mu1,mu2,mu3,z1,z2,z3,U1,U2,U3,d)
    return (
        -9/2 * mu2 * (-z2 + z1) * mu1 / (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3 + 4.5*mu3*mu1*d^3) * U1
         - 9/2 * mu2 * (2 * mu1 * z1 - 2 * mu1 * z2 - mu2 * z1) / ((-z2 + z1) * (4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) + 4.5 *mu3*mu2* d^4) * (z2 - z1) * U2
          + 3/2 * (4 * mu1 * z1 - 4 * mu1 * z2 - 3 * mu2 * z1) * mu2 / ((4 * mu1 * mu2 * z1 * z2 ^ 2 - 8 * mu1 * mu2 * z1 * z2 * z3 + 4 * mu1 * mu2 * z1 * z3 ^ 2 - 4 * mu1 * mu2 * z2 ^ 3 + 8 * mu1 * mu2 * z2 ^ 2 * z3 - 4 * mu1 * mu2 * z2 * z3 ^ 2 + 3 * mu1 * mu3 * z1 ^ 2 * z2 - 3 * mu1 * mu3 * z1 ^ 2 * z3 - 6 * mu1 * mu3 * z1 * z2 ^ 2 + 6 * mu1 * mu3 * z1 * z2 * z3 + 3 * mu1 * mu3 * z2 ^ 3 - 3 * mu1 * mu3 * z2 ^ 2 * z3 - 3 * mu2 ^ 2 * z1 * z2 ^ 2 + 6 * mu2 ^ 2 * z1 * z2 * z3 - 3 * mu2 ^ 2 * z1 * z3 ^ 2 - 3 * mu2 * mu3 * z1 ^ 2 * z2 + 3 * mu2 * mu3 * z1 ^ 2 * z3 + 3 * mu2 * mu3 * z1 * z2 ^ 2 - 3 * mu2 * mu3 * z1 * z2 * z3) * (-z3 + z2) + 4.5 *mu3*mu3* d^4) * (z3 - z2) * U3
    )
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
    circshift!(fip, state.pressure, (1, 0))
    circshift!(fjp, state.pressure, (0, 1))
    circshift!(fim, state.pressure, (-1, 0))
    circshift!(fjm, state.pressure, (0, -1))
    # Diagonal elements  
    circshift!(fipjp, state.pressure, (1, 1))
    circshift!(fimjp, state.pressure, (-1, 1))
    circshift!(fimjm, state.pressure, (-1, -1))
    circshift!(fipjm, state.pressure, (1, -1))
    # In the end it is just a weighted sum...
    @. state.h∇px .=
        state.height * (-1 / 3 * (fip - fim) - 1 / 12 * (fipjp - fimjp - fimjm + fipjm))
    @. state.h∇py .=
        state.height * (-1 / 3 * (fjp - fjm) - 1 / 12 * (fipjp + fimjp - fimjm - fipjm))

    return nothing
end

function h∇p!(state::MultiLayer_2D,sys::SysConstMultiLayer)
    fip, fjp, fim, fjm, fipjp, fimjp, fimjm, fipjm = viewneighborsMultiLayer(state.dgrad)
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
    # if sys.layers==2
    state.h∇px[:,:,1] .= state.height[:,:,1] .* (
        -1/3 .* (fip[:,:,1] .- fim[:,:,1]) .- 1/12 .* (fipjp[:,:,1] .- fimjp[:,:,1] .- fimjm[:,:,1] .+ fipjm[:,:,1])
        -1/3 .* (fip[:,:,2] .- fim[:,:,2]) .- 1/12 .* (fipjp[:,:,2] .- fimjp[:,:,2] .- fimjm[:,:,2] .+ fipjm[:,:,2])
        )
    state.h∇py[:,:,1] .= state.height[:,:,1] .* (
        -1/3 .* (fjp[:,:,1] .- fjm[:,:,1]) .- 1/12 .* (fipjp[:,:,1] .+ fimjp[:,:,1] .- fimjm[:,:,1] .- fipjm[:,:,1])
        -1/3 .* (fjp[:,:,2] .- fjm[:,:,2]) .- 1/12 .* (fipjp[:,:,2] .+ fimjp[:,:,2] .- fimjm[:,:,2] .- fipjm[:,:,2])
        )
    state.h∇px[:,:,2] .= state.height[:,:,2] .* (-1/3 .* (fip[:,:,2] .- fim[:,:,2]) .- 1/12 .* (fipjp[:,:,2] .- fimjp[:,:,2] .- fimjm[:,:,2] .+ fipjm[:,:,2]))
    state.h∇py[:,:,2] .= state.height[:,:,2] .* (-1/3 .* (fjp[:,:,2] .- fjm[:,:,2]) .- 1/12 .* (fipjp[:,:,2] .+ fimjp[:,:,2] .- fimjm[:,:,2] .- fipjm[:,:,2]))
    # end
    return nothing
end


function h∇p!(state::StateMultiLayer_1D,sys::SysConstMultiLayer_1D)
    fip, fim = viewneighborsMultiLayer_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    # In the end it is just a weighted sum...
    if sys.layers==2
        # state.h∇p[:,1] .= state.height[:,1] .* -0.5 .* (fip[:,1] .- fim[:,1] .+ fip[:,2] .- fim[:,2])
        state.h∇p[:,1] .= state.height[:,1] .* -0.5 .* (fip[:,1] .- fim[:,1] )
        state.h∇p[:,2] .= state.height[:,2] .* -0.5 .* (fip[:,2] .- fim[:,2])
    elseif sys.layers==3
        state.h∇p[:,1] .= state.height[:,1] .* -0.5 .* (fip[:,1] .- fim[:,1])
        state.h∇p[:,2] .= state.height[:,2] .* -0.5 .* (fip[:,2] .- fim[:,2])
        state.h∇p[:,3] .= state.height[:,3] .* -0.5 .* (fip[:,3] .- fim[:,3])
    end
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

# for Miscible films
function h∇p!(state::StateMiscible_1D,sys::SysConstMiscible_1D)
    fip, fim = viewneighborsMiscible_1D(state.dgrad)
    # One dim case, central differences
    circshift!(fip, state.pressure, 1)
    circshift!(fim, state.pressure, -1)
    # In the end it is just a weighted sum...
    # if sys.layers==2
        state.h∇p .= state.height .* -0.5 .* (fip .- fim)
    # end
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
    circshift!(fip, state.basestate.pressure, (1, 0))
    circshift!(fjp, state.basestate.pressure, (0, 1))
    circshift!(fim, state.basestate.pressure, (-1, 0))
    circshift!(fjm, state.basestate.pressure, (0, -1))
    # Diagonal elements  
    circshift!(fipjp, state.basestate.pressure, (1, 1))
    circshift!(fimjp, state.basestate.pressure, (-1, 1))
    circshift!(fimjm, state.basestate.pressure, (-1, -1))
    circshift!(fipjm, state.basestate.pressure, (1, -1))
    # In the end it is just a weighted sum...
    @. state.basestate.h∇px .=
        state.basestate.height *
        (-1 / 3 * (fip - fim) - 1 / 12 * (fipjp - fimjp - fimjm + fipjm))
    @. state.basestate.h∇py .=
        state.basestate.height *
        (-1 / 3 * (fjp - fjm) - 1 / 12 * (fipjp + fimjp - fimjm - fipjm))

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
    ∇gamma!(state)

Computes the gradient of a spatially resolved surface tension field via central differneces and then concludes the surface stress as 

`` \\frac{h^2 + 2 \\delta h }{2M(h)}``

with 

``M(h)= \\frac{2h^2 + 6 \\delta h + 3\\delta^2 }{ 6}``

and implicitly ``\\rho_h=1``. 
"""
function ∇gamma!(state::StateMiscible_1D, sys::SysConstMiscible_1D)
    fip, fim = viewneighbors_1D(state.dgrad[:,:,1])
    # One dim case, central differences
    circshift!(fip, state.gamma[:,1], 1)
    circshift!(fim, state.gamma[:,1], -1)
    
    # In the end it is just a weighted sum...
    # state.∇γ .= -3/2 .* ((fip .- fim) ./ 2.0)
    state.∇gamma .=((fip .- fim) .* -0.5)
    state.stress[:,1] .= state.∇gamma .* ((state.height[:,1] .+ state.height[:,2]) .* (state.height[:,1] .+ state.height[:,2]) .* 0.5 .+ sys.delta .*(state.height[:,1] .+ state.height[:,2])) .* 6 ./ (2 .* (state.height[:,1] .+ state.height[:,2]) .* (state.height[:,1] .+ state.height[:,2]) .+ 6*sys.delta .* (state.height[:,1] .+ state.height[:,2]) .+ 3*sys.delta*sys.delta )
    state.stress[:,2] .= state.stress[:,1]
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
    fluc_x .*=
        sqrt.(
            2 .* kᵦT .* μ .* 6 .* height ./
            (2 .* height .* height .+ 6 .* height .* δ .+ 3 .* δ .* δ)
        )
    fluc_y .*=
        sqrt.(
            2 .* kᵦT .* μ .* 6 .* height ./
            (2 .* height .* height .+ 6 .* height .* δ .+ 3 .* δ .* δ)
        )
    return nothing
end

thermal!(state::State_thermal, sys::SysConst) = thermal!(
    state.kbtx,
    state.kbty,
    state.basestate.height,
    sys.param.kbt,
    sys.param.μ,
    sys.param.δ,
)

function thermal!(fluc, height, kᵦT, μ, δ)
    randn!(fluc)
    fluc .*=
        sqrt.(
            2 .* kᵦT .* μ .* 6 .* height ./
            (2 .* height .* height .+ 6 .* height .* δ .+ 3 .* δ .* δ)
        )

    return nothing
end

thermal!(state::State_thermal_1D, sys::Consts_1D) =
    thermal!(state.kbt, state.basestate.height, sys.param.kbt, sys.param.μ, sys.param.δ)

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

See also: [`Swalbe.run_dropletforced`](@ref)
"""
function inclination!(α::Vector, state::LBM_state_2D; t = 1000, tstart = 0, tsmooth = 1)
    @. state.Fx .+= state.height * α[1] * (0.5 + 0.5 * tanh((t - tstart) / tsmooth))
    @. state.Fy .+= state.height * α[2] * (0.5 + 0.5 * tanh((t - tstart) / tsmooth))

    return nothing
end

function inclination!(α::Vector, state::Expanded_2D; t = 1000, tstart = 0, tsmooth = 1)
    @. state.basestate.Fx .+=
        state.basestate.height * α[1] * (0.5 + 0.5 * tanh((t - tstart) / tsmooth))
    @. state.basestate.Fy .+=
        state.basestate.height * α[2] * (0.5 + 0.5 * tanh((t - tstart) / tsmooth))

    return nothing
end

function inclination!(α::Float64, state::State_1D; t = 1000, tstart = 0, tsmooth = 1)
    state.F .+= state.height .* α .* (0.5 .+ 0.5 .* tanh((t - tstart) / tsmooth))

    return nothing
end

function inclination!(α::Float64, state::Expanded_1D; t = 1000, tstart = 0, tsmooth = 1)
    state.basestate.F .+=
        state.basestate.height .* α .* (0.5 .+ 0.5 .* tanh((t - tstart) / tsmooth))

    return nothing
end

"""
    update_rho()

Time evolution of the `active` field rho.

TODO: @Tilman!
"""
function update_rho!(rho, rho_int, height, dgrad, differentials; D = 1.0, M = 0.0)
    lap_rho, grad_rho, lap_h, grad_h = view_four(differentials)
    ∇²f!(lap_rho, rho, dgrad)
    ∇f!(grad_rho, rho, dgrad)
    ∇²f!(lap_h, height, dgrad)
    ∇f!(grad_h, height, dgrad)
    rho_int .= (
        D .* lap_rho .- M .* (grad_rho .^ 2 .+ rho .* lap_rho) .-
        D .* (
            grad_rho .* grad_h ./ height +
            rho .* (lap_h ./ height .- (grad_h ./ height) .^ 2)
        )
    )

    rho .+= rho_int

    return nothing
end

"""
    surface_tension_gradient!(state)

Computes the gradient of a spatially resolved surface tension field.
"""
function ∇γ!(state::T) where {T<:Expanded_1D}
    fip, fim = viewneighbors_1D(state.basestate.dgrad)
    # One dim case, central differences
    circshift!(fip, state.γ, 1)
    circshift!(fim, state.γ, -1)

    # In the end it is just a weighted sum...
    state.∇γ .= -3 / 2 .* ((fip .- fim) ./ 2.0)
    return nothing
end

function ∇γ!(state::T, sys::SysConst_1D) where {T<:Expanded_1D}
    fip, fim = viewneighbors_1D(state.basestate.dgrad)
    # One dim case, central differences
    circshift!(fip, state.γ, 1)
    circshift!(fim, state.γ, -1)

    # Here I try to get rid of the slippage term in the gradient
    state.∇γ .=
        (
            2 .* state.basestate.height .^ 2 .+
            6 .* sys.param.δ .* state.basestate.height .+ 3 * sys.param.δ^2
        ) ./ (6 .* state.basestate.height) .* state.basestate.height ./ 2 .*
        ((fip .- fim) ./ 2.0)
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

