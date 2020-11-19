@testset "Collision and Stream" begin

end

@testset "viewdists" begin
    f = reshape(collect(1.0:225.0),5,5,9)
    f0, f1, f2, f3, f4, f5, f6, f7, f8 = Swalbe.viewdists(f)
    allviews = [f0, f1, f2, f3, f4, f5, f6, f7, f8]
    for (index, value) in enumerate(allviews)
        @test all(value .== f[:,:,index])
    end
end