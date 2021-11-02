# Swalbe.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://zitzeronion.github.io/Swalbe.jl/dev/)
[![CI](https://github.com/Zitzeronion/Swalbe.jl/workflows/CI/badge.svg?branch=master&event=push)](https://github.com/Zitzeronion/Swalbe.jl/actions)
[![codecov](https://codecov.io/gh/Zitzeronion/Swalbe.jl/branch/master/graph/badge.svg?token=J1AMK7YW69)](https://codecov.io/gh/Zitzeronion/Swalbe.jl)
[![status](https://joss.theoj.org/papers/414a5b53a41e05a250a352360a7da337/status.svg)](https://joss.theoj.org/papers/414a5b53a41e05a250a352360a7da337)

![Dewetting_logo](https://gist.githubusercontent.com/Zitzeronion/807b9a7b2226e65643288df9a8cc1f46/raw/3a561e2a2b09eb42bf688f1d304f658b93fba8ed/logo_animation.gif)


## Thin film simulations using lattice Boltzmann :rainbow: :ocean:

Why is a thin film solver called **Swalbe.jl** you may ask?

The idea is to use the
[*lattice Boltzmann method (LBM)*](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods)
and all its benefits (easy to code, vast amount of literature and scalability)
to simulate thin liquid film flows. Instead of reinventing the wheel we make use 
of a class of lattice Boltzmann models that were build to simulate shallow water 
problems, see
  [Salmon](http://pordlabs.ucsd.edu/rsalmon/salmon.1999a.pdf) (not the fish :fish:),
  [Dellar](https://people.maths.ox.ac.uk/dellar/papers/LBshallow.pdf) and
  [van Thang et al.](https://hal.archives-ouvertes.fr/hal-01625073/document) (*all free to read*).
Thus the name of the package **S**hallow **WA**ter **L**attice **B**oltzmann slov**E**r or **Swalbe**.

Of course using a plain shallow water model will not work to simulate thin film
dynamics, that is the reason we build our own model :neckbeard:.  Now the main
difference is that we throw away most of the shallow water parts by assuming
they are small as compared to thin film relevant things, e.g. the substrate
fluid interaction.  The full explanation of the model with some benchmarks can
be found in our paper
[Zitz et al.](http://pub.hi-ern.de/publications/2019/ZSMDH19/2019-ThinFilm-PRE.pdf)
(the C/C++ OpenACC codebase has not been further developed since the project
moved to *Julia*)

## How to **get**

### Requirements
First of all you need a *Julia* (>= 1.6) installation.  *Julia* is a high level
open source programming language and it is as easy to use as python :snake: (my
opinion).

*Julia* can be downloaded at the projects homepage
[julialang.org](https://julialang.org/).

### Install using the Julia package manager

**Swalbe.jl** is a registered package of the *Julia* package manager.  The only
thing you have to do is to add the package to your *Julia* environment with:

```julia
julia> ] add Swalbe
```

### Install from source

Of course you can also clone or fork the repo and activate the package inside
the julia **REPL** (Read Evaluate Print Loop). First you need to go the Swalbe
directory and open a **REPL**

```bash
git clone <swalbe git url>
cd swalbe.jl
julia
```

now you can activate the package with

```julia
julia> ] activate .
```

To check if the package works you can run the test suite with

```julia
julia> ] test Swalbe
```

All tests can be found in
[test folder](https://github.com/Zitzeronion/Swalbe.jl/tree/master/test), but do not
expect too many comments. Still especially the
[simulate.jl](https://github.com/Zitzeronion/Swalbe.jl/blob/master/test/simulate.jl)
file is worth a look.

## How to **use**

The idea of **Swalbe.jl** is to script your thin film simulation, based on a lattice Boltzmann iteration.  
That is why most core functions can be easily extended, or used out of the box. Find some examples in the [scripts](https://github.com/Zitzeronion/Swalbe.jl/tree/master/scripts) folder.

Some initial conditions are handily pre-programmed. 
E.g. simulating the Rayleigh-Taylor instability:


```julia
using Swalbe

# set the constants of the system
sys = Swalbe.SysConst(Lx=100, Ly=100, g=-0.001, γ=0.0005, Tmax=1000)

# run with given parameters for Tmax timesteps ...
# return Lx×Ly array with the final configuration
h = Swalbe.run_rayleightaylor(sys, "CPU"; h₀=1.0, ϵ=0.01, verbos=true)
```

Further examples can be found in the tutorials section of the documentation: [Tutorials](https://zitzeronion.github.io/Swalbe.jl/dev/tutorials/)

## How to **support and contribute**

Leave a star if you like the idea of the project and/or the content of the
package.  You can support the project by actively using it and raising
[issues](https://github.com/Zitzeronion/Swalbe.jl/issues).
Help is always very welcome, if you want to contribute open a
[PR](https://github.com/Zitzeronion/Swalbe.jl/pulls) or raise an
[issue](https://github.com/Zitzeronion/Swalbe.jl/issues) with a feature request
(and if possible with a way how to include it).  Feel free to DM me on
[Twitter](https://twitter.com/Zitzero) if you have questions, I will try to
answer them all timely.

## **Status** of the Package

![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)

<!-- ![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->

The package has reached a stable release with version 0.1.
All tools needed for running a numerical experiment are tested and usable.
I am currently writing a paper for which all experiments were done with this package :blush:.

## **Credit**

- [Stefan Zitz (main developer)](https://www.hi-ern.de/hi-ern/CompFlu/Team/Zitz/zitz.html)
- [Manuel Zellhöfer](https://www.hi-ern.de/hi-ern/CompFlu/Team/Zellhoefer/zellhoefer.html)
- [Andrea Scagliarini](https://www.iac.rm.cnr.it/iacsite/index.php?page=people&id=140)
- [Jens Harting](https://www.hi-ern.de/hi-ern/CompFlu/Team/Harting/harting.html)