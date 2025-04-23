# Swalbe.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://zitzeronion.github.io/Swalbe.jl/stable/)
[![CI](https://github.com/Zitzeronion/Swalbe.jl/actions/workflows/ci.yml/badge.svg?branch=master&event=push)](https://github.com/Zitzeronion/Swalbe.jl/actions)
[![codecov](https://codecov.io/gh/Zitzeronion/Swalbe.jl/branch/master/graph/badge.svg?token=J1AMK7YW69)](https://codecov.io/gh/Zitzeronion/Swalbe.jl)
[![status](https://joss.theoj.org/papers/414a5b53a41e05a250a352360a7da337/status.svg)](https://joss.theoj.org/papers/414a5b53a41e05a250a352360a7da337)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7030890.svg)](https://doi.org/10.5281/zenodo.7030890)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FSwalbe&query=total_requests&label=Total)](http://juliapkgstats.com/pkg/Swalbe) 
<!--[![Swalbe.jl Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Swalbe)](https://pkgs.genieframework.com?packages=Swalbe) --->




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

Some development for this solver was performed under the priority program **[SPP2171-Dynamic Wetting of Flexible, Adaptive, and Switchable Surfaces](https://www.uni-muenster.de/SPP2171/index.html)**. 
On the homepage of the SPP in the resources' section we supply a simple tutorial for the [coalescence of droplets](https://www.uni-muenster.de/imperia/md/content/SPP2171/droplet_coalescence_tutorial.html) using a [**Pluto**](https://github.com/fonsp/Pluto.jl) notebook.

## How to **support and contribute**

Leave a star if you like the idea of the project and/or the content of the
package.  You can support the project by actively using it and raising
[issues](https://github.com/Zitzeronion/Swalbe.jl/issues).
Help is always very welcome, if you want to contribute open a
[PR](https://github.com/Zitzeronion/Swalbe.jl/pulls) or raise an
[issue](https://github.com/Zitzeronion/Swalbe.jl/issues) with a feature request
(and if possible with a way how to include it).  Feel free to DM me on
[Bluesky](https://bsky.app/profile/zitzeronion.bsky.social) if you have questions, I will try to
answer them all timely.

## Publications

As of April 2025 this package has been used to generate 6 peer reviewed publications which can be found in the [docs](https://zitzeronion.github.io/Swalbe.jl/stable/).
Furthermore there is a published doctoral [thesis](https://open.fau.de/items/bc35f5d4-b157-4aad-82be-9d23825cb37a/full) and one that is soon to be submitted tied to this projects, link supplied in the docs. 

## **Credit**

If you use this package please cite it as
```bibtex
@article{Zitz2022,
  doi = {10.21105/joss.04312},
  url = {https://doi.org/10.21105/joss.04312},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {77},
  pages = {4312},
  author = {Stefan Zitz and Manuel Zellhöfer and Andrea Scagliarini and Jens Harting},
  title = {Swalbe.jl: A lattice Boltzmann solver for thin film hydrodynamics},
  journal = {Journal of Open Source Software}
} 
```

- [Stefan Zitz (main developer)](https://forskning.ruc.dk/en/persons/zitz)
- [Tilman Richter](https://www.hi-ern.de/profile/richter_t)
- [Manuel Zellhöfer](https://www.hi-ern.de/hi-ern/CompFlu/Team/Zellhoefer/zellhoefer.html)
- [Andrea Scagliarini](https://www.iac.rm.cnr.it/iacsite/index.php?page=people&id=140)
- [Jens Harting](https://www.hi-ern.de/de/forschung/dynamik-komplexer-fluide-und-grenzflaechen)
