module Swalbe
# Coding
import Revise, CUDA, Statistics, Random, Base, LazySets 
# Datatypes and output
import FileIO, JLD2, BSON, DataFrames

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

end
