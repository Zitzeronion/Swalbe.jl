using Swalbe, Test, Statistics, Parameters, Random, LazySets

@testset "Swalbe.jl" begin
    include("initialize.jl")
    include("initialvalues.jl")
    include("differences.jl")
    include("pressure.jl")
    include("collide.jl")
    include("equilibrium.jl")
    include("forcing.jl")
    include("moments.jl")
end
