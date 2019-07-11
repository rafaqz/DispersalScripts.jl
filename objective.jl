struct OffsetRegionObjective{DT,RL,OC,FS,S} <: AbstractObjective
    detectionthreshold::DT
    regionlookup::RL
    occurance::OC
    framesperstep::FS
    start::S
end

Dispersal.targets(obj::OffsetRegionObjective) = obj.occurance

Dispersal.simpredictions(obj::OffsetRegionObjective, simoutput) = begin
    nregions, nsteps = size(obj.occurance)
    nframes = length(simoutput)
    nsteps = ceil(Int, nframes / obj.framesperstep)
    # Allocate arrays for steps and set all cells to zero 
    steps = [zeros(Bool, size(simoutput[1])) for f in 1:nsteps]
    fill!.(steps, zero(eltype(steps[1])))

    # Get the mean population for steps from the frames in each step:
    # First incomplete step
    firststep = 1:obj.framesperstep - obj.start + 1
    for frame in firststep  
        steps[1] .|= poptoobs.(simoutput[frame], obj.detectionthreshold)
    end
    # Full steps
    remainingsteps = last(firststep) + 1:nframes
    frame = 3
    for frame in remainingsteps
        step = stepfromframe(obj.framesperstep, frame, obj.start)
        steps[step] .|= poptoobs.(simoutput[frame], obj.detectionthreshold)
    end

    # Allocate a boolean array to contain our presence/absence predictions
    prediction = zeros(Bool, size(obj.occurance))

    # Convert mean cell populations to regional prescence/absence
    for s in 1:nsteps
        for r in 1:nregions
            reg = celltoregion(obj.regionlookup, steps[s], r)
            prediction[r, s] = reg
        end
    end
    prediction 
end


# Which cells are above the detection threshold for observations
poptoobs(pop, threshold) = pop >= threshold

# Do any cells in the region contain observations
celltoregion(regionlookup::AbstractArray, cells::AbstractArray, region_id) = 
    any((regionlookup .== region_id) .& cells)

stepfromframe(framesperstep, start, t) = (t - 2one(t) + start) ÷ framesperstep + one(t)


# Tests

using Test

@testset "Offset Region Objective" begin

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

    @test poptoobs.(pop, 2) == Bool[0 1 1
                                    0 0 1
                                    1 1 1]

    @test poptoobs.(pop, 3) == Bool[0 0 1
                                    0 0 1
                                    0 0 1]

    # Test celltoregion() for five regions using 2 thresholds
    lookup = [2 2 3 
              1 3 3
              4 5 5]

    step = poptoobs.(pop, 2)
    @test celltoregion.(Ref(lookup), Ref(step), 1:5) == [false, true, true, true, true]

    step = poptoobs.(pop, 3)
    @test celltoregion.(Ref(lookup), Ref(step), 1:5) == [false, false, true, false, true]

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
    objective = OffsetRegionObjective(detectionthreshold, lookup, occurance, framesperstep, start)
    obj = objective
    simoutput = output
    prediction = Bool[0 0 0 
                      0 0 1
                      0 1 0
                      0 1 1
                      1 1 1]

    @test prediction == Dispersal.simpredictions(objective, output)
end
