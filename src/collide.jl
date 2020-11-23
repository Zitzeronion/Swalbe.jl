"""
    BGKandStream!(fout, feq, ftemp, Fx, Fy, τ)

Performs a BGK collision operation with a WFM forcecorrection and a subsequent streaming of the resulting populations.

# Arguments

- `fout`: Streamed distribution after the collision processes.
- `feq`: Equilibrium distribution, computed with equilibrium!.
- `ftemp`: Temporary distribution from the time step before, only useful if `` \\tau \\neq 1``.
- `Fx`: Sum of forces acting on the fluid in x-direction.
- `Fy`: Sum of forces acting on the fluid in y-direction.
- `τ`: Relaxtion time.

# Mathematics

The lattice Boltzmann equation in its discretized format is relatively simple to write down

`` f_{\\alpha}(\\mathbf{x}+\\mathbf{e}_{\\alpha}\\Delta t, t+\\Delta t) - f_{\\alpha}(\\mathbf{x}, t) = -\\frac{1}{\\tau}(f_{\\alpha}(\\mathbf{x}, t) - f^{\\text{eq}}_{\\alpha}(\\mathbf{x}, t)) + \\Delta t \\mathcal{S}_{\\alpha}, ``

where the collision kernel is approximated with a BKG single relaxation time (SRT) 

`` \\Omega_{\\alpha} = \\frac{1}{\\tau}(f_{\\alpha}(\\mathbf{x}, t) - f^{\\text{eq}}_{\\alpha}(\\mathbf{x}, t), ``

and a source term `` \\mathcal{S} `` which is given by

`` \\mathcal{S}_{\\alpha} = \\frac{3 w_{\\alpha}}{e_{\\alpha x}e_{\\alpha x}+e_{\\alpha y}e_{\\alpha y}}\\mathbf{e_{\\alpha}}\\cdot\\mathbf{F}_{\\alpha} . ``

The term `` e_{\\alpha x}e_{\\alpha x}+e_{\\alpha y}e_{\\alpha y} `` is either zero for the zeroth population, 1 for the first four populations or 2 for the remaining ones.

# Examples
    
```jldoctest
julia> using Swalbe

julia> feq = ones(5,5,9); ftemp = zeros(5,5,9); fout = zeros(5,5,9);

julia> feq[1,1,:] .= 2.0 # To check the streaming process 
9-element view(::Array{Float64,3}, 1, 1, :) with eltype Float64:
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0

julia> Swalbe.BGKandStream!(fout, feq, ftemp, zeros(5,5), zeros(5,5))

julia> fout[:,:,6] # The value 2 should have moved one down and one to the right!
5×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
 1.0  2.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0

```

# References

- [Salmon](https://www.ingentaconnect.com/contentone/jmr/jmr/1999/00000057/00000003/art00005#)
- [Dellar](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.036309)
- [Peng et al.](https://onlinelibrary.wiley.com/doi/full/10.1002/fld.4726)

See also: [`Swalbe.equilibrium`](@ref)
"""
function BGKandStream!(fout, feq, ftemp, Fx, Fy, τ)
    # All distribution functions
    fe0, fe1, fe2, fe3, fe4, fe5, fe6, fe7, fe8 = viewdists(feq)
    ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8 = viewdists(ftemp)
    fo0, fo1, fo2, fo3, fo4, fo5, fo6, fo7, fo8 = viewdists(fout)

    omeg = 1-1/τ
    # Collision for the nine populations
    # Zeroth, no forcing!
    fo0 .= omeg .* ft0 .+ 1/τ .* fe0
    # Straight ones, force correction is 1/3, 3*1/9
    fo1 .= omeg .* ft1 .+ 1/τ .* fe1 .+ 1/3 .* Fx
    fo2 .= omeg .* ft2 .+ 1/τ .* fe2 .+ 1/3 .* Fy
    fo3 .= omeg .* ft3 .+ 1/τ .* fe3 .- 1/3 .* Fx
    fo4 .= omeg .* ft4 .+ 1/τ .* fe4 .- 1/3 .* Fy
    # Diagonal ones, force correction 1/24 -> 3/2*1/36
    fo5 .= omeg .* ft5 .+ 1/τ .* fe5 .+ 1/24 .* (Fx .+ Fy)
    fo6 .= omeg .* ft6 .+ 1/τ .* fe6 .+ 1/24 .* (Fy .- Fx)
    fo7 .= omeg .* ft7 .+ 1/τ .* fe7 .- 1/24 .* (Fx .+ Fy)
    fo8 .= omeg .* ft8 .+ 1/τ .* fe8 .+ 1/24 .* (Fx .- Fy)

    # This is the streaming step with implicite periodic boundarys
    ft0 .= circshift(fo0, (0, 0))
    ft1 .= circshift(fo1, (1, 0))
    ft2 .= circshift(fo2, (0, 1))
    ft3 .= circshift(fo3, (-1, 0))
    ft4 .= circshift(fo4, (0, -1))
    ft5 .= circshift(fo5, (1, 1))
    ft6 .= circshift(fo6, (-1, 1))
    ft7 .= circshift(fo7, (-1, -1))
    ft8 .= circshift(fo8, (1, -1))
    
    # Overwrite fout with ftemp
    fout .= ftemp
    return nothing
end

# Explicit version for the case τ = 1, which simplifies the computation quite a lot.
function BGKandStream!(fout, feq, ftemp, Fx, Fy)
    # All distribution functions
    fe0, fe1, fe2, fe3, fe4, fe5, fe6, fe7, fe8 = viewdists(feq)
    ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8 = viewdists(ftemp)
    fo0, fo1, fo2, fo3, fo4, fo5, fo6, fo7, fo8 = viewdists(fout)

    # Collision for the nine populations
    # Zeroth, no forcing!
    fo0 .= fe0
    # Straight ones, force correction is 1/3, 3*1/9
    fo1 .= fe1 .+ 1/3 .* Fx
    fo2 .= fe2 .+ 1/3 .* Fy
    fo3 .= fe3 .- 1/3 .* Fx
    fo4 .= fe4 .- 1/3 .* Fy
    # Diagonal ones, force correction 1/24 -> 3/2*1/36
    fo5 .= fe5 .+ 1/24 .* (Fx .+ Fy)
    fo6 .= fe6 .+ 1/24 .* (Fy .- Fx)
    fo7 .= fe7 .- 1/24 .* (Fx .+ Fy)
    fo8 .= fe8 .+ 1/24 .* (Fx .- Fy)

    # This is the streaming step with implicite periodic boundarys
    ft0 .= circshift(fo0, (0, 0))
    ft1 .= circshift(fo1, (1, 0))
    ft2 .= circshift(fo2, (0, 1))
    ft3 .= circshift(fo3, (-1, 0))
    ft4 .= circshift(fo4, (0, -1))
    ft5 .= circshift(fo5, (1, 1))
    ft6 .= circshift(fo6, (-1, 1))
    ft7 .= circshift(fo7, (-1, -1))
    ft8 .= circshift(fo8, (1, -1))
    
    # Overwrite fout with ftemp
    fout .= ftemp
    return nothing
end

"""
    viewdists(f)

Generates a view for all nine populations of a **D2Q9** distribution function.

# Examples
```jldoctest
julia> using Swalbe, Test

julia> ftest = reshape(collect(1.0:225.0),5,5,9);

julia> f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(ftest);

julia> @test all(f3 .== ftest[:,:,4])
Test Passed

```

See also: [`Swalbe.BGKandStream`](@ref)
"""
function viewdists(f)
    f0 = view(f, :, :, 1)
    f1 = view(f, :, :, 2)
    f2 = view(f, :, :, 3)
    f3 = view(f, :, :, 4)
    f4 = view(f, :, :, 5)
    f5 = view(f, :, :, 6)
    f6 = view(f, :, :, 7)
    f7 = view(f, :, :, 8)
    f8 = view(f, :, :, 9)
    
    return f0, f1, f2, f3, f4, f5, f6, f7, f8
end