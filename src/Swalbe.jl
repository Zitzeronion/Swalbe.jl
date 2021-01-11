module Swalbe

using Revise, CUDA, Statistics, Parameters, Random, LazySets

include("initialize.jl")
include("initialvalues.jl")
include("differences.jl")
include("pressure.jl")
include("collide.jl")
include("equilibrium.jl")
include("forcing.jl")
include("moments.jl")
include("simulate.jl")
include("experiments.jl")

end
