using Swalbe, Test, LinearAlgebra, DiffEqOperators, Images

@testset "Swalbe.jl" begin
    include("differences.jl")
    include("pressure.jl")
end
