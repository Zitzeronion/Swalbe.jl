# Swalbe.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://zitzeronion.github.io/Swalbe.jl/dev/) [![CI](https://github.com/Zitzeronion/Swalbe.jl/workflows/CI/badge.svg?branch=master&event=push)](https://github.com/Zitzeronion/Swalbe.jl/actions) [![codecov](https://codecov.io/gh/Zitzeronion/Swalbe.jl/branch/master/graph/badge.svg?token=J1AMK7YW69)](https://codecov.io/gh/Zitzeronion/Swalbe.jl)

## Thin film simulations using lattice Boltzmann :rainbow: :ocean:

Why is a thin film solver called **Swalbe.jl** you may ask?

The idea is to use the [*lattice Boltzmann method (LBM)*](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) to simulate thin liquid film flows.
Instead of reinventing the wheel we make use of a class of lattice Boltzmann models that were build to simulate shallow water problems, see [Salmon](http://pordlabs.ucsd.edu/rsalmon/salmon.1999a.pdf)(not the fish :fish:), [Dellar](https://people.maths.ox.ac.uk/dellar/papers/LBshallow.pdf) and [van Thang et al.](https://hal.archives-ouvertes.fr/hal-01625073/document) (*all free to read*).
Thus the name of the package **S**hallow **WA**ter **L**attice **B**oltzmann slov**E**r or **Swalbe**.

Of course using a plain shallow water model will not work to simulate thin film dynamics, that is the reason we build our own model :neckbeard:

## How to **get** 

First of all you need a *Julia* installation. 
*Julia* is a high level open source programming language and it is as easy to use as python :snake: (my opinion).

*Julia* can be downloaded at the projects homepage [julialang.org](https://julialang.org/), use the correct installation for your operating system.
**Important** for **CUDA** we require *Julia* version 1.6 or higher.  

Since recently **Swalbe.jl** is a registered package of the *Julia* package manager.
The only thing you have to do is to add the package to your *Julia* environment with: 
```julia
julia> ] add Swalbe
```

To see that the package works you can run the test suit with
```julia
julia> ] test Swalbe
```
All tests can be found in test folder, but do not expect too many comments there.