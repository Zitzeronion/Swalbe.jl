using Swalbe, Test, Statistics, Parameters, Random

@testset "Swalbe.jl" begin
    include("initialize.jl")
    include("differences.jl")
    include("pressure.jl")
    include("collide.jl")
    include("equilibrium.jl")
    include("forcing.jl")
end
