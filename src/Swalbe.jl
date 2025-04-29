module Swalbe
# Coding
using Revise, CUDA, Statistics, Random, Base, LazySets
# Datatypes and output
import FileIO, BSON 
import DataFrames: DataFrame
import JLD2: load, save


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
