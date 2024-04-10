module Swalbe
# Coding
using Revise, Statistics, Random, Base, LazySets, ThreadsX, LinearAlgebra, FFTW
using CUDA
# Datatypes and output
using FileIO, JLD2, BSON, DataFrames

include("initialize.jl")
include("initialvalues.jl")
include("differences.jl")
include("pressure.jl")
include("collide.jl")
include("equilibrium.jl")
include("forcing.jl")
include("moments.jl")
include("simulate.jl")
include("measures.jl")
include("logic.jl")
include("active.jl")

end
