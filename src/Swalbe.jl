module Swalbe

using Revise, CUDA, Statistics, Parameters, Random

include("initialize.jl")
include("differences.jl")
include("pressure.jl")
include("collide.jl")
include("equilibrium.jl")
include("forcing.jl")

end
