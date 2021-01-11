@testset "Measurements" begin
    @testset "Wetted liquid solid" begin
        area_size = [1.0]
        maxheight = [1.0]
        height = zeros(5,5)
        height[2:3, 2:3] .= 1.0
        Swalbe.wetted!(area_size, maxheight, height, 1)
        @test area_size[1] == 4.0
        @test maxheight[1] == 1.0
    
        area_size = ones(5)
        height = zeros(5,5)
        maxheight = zeros(5)
        for t in 1:5
            height[t, t] = t
            Swalbe.wetted!(area_size, maxheight, height, t)
        end
        sol1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        for t in 1:5
            @test area_size[t] == sol1[t]
            @test maxheight[t] == sol1[t]
        end
    end
    @testset "Position liquid" begin
        height = zeros(5,5)
        fluid = falses(1,25)
        height[2:3, 2:3] .= 1.0
        dummy = falses(5,5)
        Swalbe.fluid_dry!(fluid, dummy, height, 1)

        @test all(fluid[1, 1:6] .== false)
        @test all(fluid[1, 7:8] .== true)
    
        fluid = falses(1,25)
        height = zeros(5,5)
        maxheight = zeros(5)
        # for t in 1:5
        #     height[t, t] = t
        #     Swalbe.wetted!(area_size, drop_pos, maxheight, height, t; hpfilm = 0.055)
        # end
        # sol1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        # for t in 1:5
        #     @test area_size[t] == sol1[t]
        #     @test maxheight[t] == sol1[t]
        # end
    end
end