using Swalbe, Test, Statistics, Parameters

@testset "Swalbe.jl" begin
    include("initialize.jl")
    include("differences.jl")
    include("pressure.jl")
    include("collide.jl")
    include("equilibrium.jl")
    include("forcing.jl")
end
