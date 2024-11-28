module Swalbe
# Coding
using Revise, CUDA, Statistics, Random, Base, LazySets, FFTW
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
include("obstacle.jl")
include("logic.jl")
include("active.jl")



end
