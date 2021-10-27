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

function moments!(state::State, device::String)
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