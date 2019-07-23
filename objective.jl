# Do any cells in the region contain observations
# celltoregion!(regionbuffer, regionlookup::AbstractArray, cells::AbstractArray, region_id) =
    # regionbuffer .= (regionlookup .== region_id) .& cells

# Tests

using Test

@testset "Region Objective" begin

    @test stepfromframe(3, 2, 1) == 1
    @test stepfromframe(3, 2, 2) == 1
    @test stepfromframe(3, 2, 3) == 2
    @test stepfromframe(3, 2, 4) == 2
    @test stepfromframe(3, 2, 5) == 2
    @test stepfromframe(3, 2, 6) == 3
    @test stepfromframe(3, 2, 7) == 3
    @test stepfromframe(3, 2, 8) == 3

    obs = zeros(Bool, 3, 3)
    pop = [1.0 2.5 5.7
           0.0 0.7 4.3
           2.3 2.9 3.0]

    # @test poptoobs.(pop, 2) == Bool[0 1 1
    #                                 0 0 1
    #                                 1 1 1]

    # @test poptoobs.(pop, 3) == Bool[0 0 1
    #                                 0 0 1
    #                                 0 0 1]

    # # Test celltoregion() for five regions using 2 thresholds
    # lookup = [2 2 3
    #           1 3 3
    #           4 5 5]

    # buf = zeros(Bool, size(lookup))

    # step = poptoobs.(pop, 2)
    # @test any.(celltoregion!.(Ref(buf), Ref(lookup), Ref(step), 1:5)) == [false, true, true, true, true]

    # step = poptoobs.(pop, 3)
    # @test any.(celltoregion!.(Ref(buf), Ref(lookup), Ref(step), 1:5)) == [false, false, true, false, true]

    # Only used for dimensions, not data
    occurance = [0 0 0
                 0 0 0
                 0 0 0
                 0 0 0
                 0 0 0]

    output = [
              [0.0 0.0 0.0
               0.0 0.0 0.0
               0.0 0.0 0.3],
              [0.0 0.0 0.0
               0.0 0.0 0.4
               0.0 0.3 2.3],
              [0.0 0.0 0.1
               0.0 0.4 0.7
               0.1 1.2 3.3],
              [0.0 0.0 0.2
               0.0 0.8 0.9
               0.7 1.8 4.0],
              [0.0 0.2 0.5
               0.0 1.0 1.0
               1.0 2.0 4.0],
              [0.0 1.0 0.2
               0.0 0.7 0.9
               2.0 3.0 4.0],
              [0.0 0.2 0.0
               0.0 0.5 0.7
               2.0 3.0 4.0],
              [0.0 0.0 0.0
               0.0 0.2 0.5
               2.0 3.0 4.0]
             ]

    lookup = [2 2 3
              1 3 3
              4 5 5]

    # Test OffsetRegionObjective
    framesperstep = 3
    detectionthreshold = 1.0
    start = 2
    objective = RegionObjective(detectionthreshold, lookup, occurance, framesperstep, start)
    obj = objective
    simoutput = output
    prediction = Bool[0 0 0
                      0 0 1
                      0 1 0
                      0 1 1
                      1 1 1]

    @test prediction == Dispersal.simpredictions(objective, output)
end
