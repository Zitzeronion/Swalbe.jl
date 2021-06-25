# Internal Functions

```@meta
CurrentModule = Swalbe
```

## Lattice Boltzmann core

```@docs
BGKandStream!(fout, feq, ftemp, Fx, Fy, τ)
```

```@docs
moments!(height, velx, vely, fout)
```


## Pressure and forces

The pressure function is at the heart of this method

```@docs
filmpressure!(pressure, height, γ, θ, n, m, hmin, hcrit)
```

- link to [Swalbe.jl Documentation](@ref)
- link to [`filmpressure!(pressure, height, γ, θ, n, m, hmin, hcrit)`](@ref)
