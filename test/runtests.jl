using Swalbe, Test, Statistics, Parameters, Random, LazySets

@testset "Swalbe.jl" begin
    println("Testing allocations")
    include("initialize.jl")
    println("Testing initial states")
    include("initialvalues.jl")
    println("Testing finite differences")
    include("differences.jl")
    println("Testing film pressure")
    include("pressure.jl")
    println("Testing collision operator")
    include("collide.jl")
    println("Testing equilibrium distribution")
    include("equilibrium.jl")
    println("Testing forces")
    include("forcing.jl")
    println("Testing moment calculation")
    include("moments.jl")
    println("Testing measurements")
    include("measures.jl")
    # println("Testing with relevant simulations")
    # include("simulate.jl")
end
