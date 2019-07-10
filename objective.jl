
struct OffsetRegionObjective{DT,RL,OC,FS,S} <: AbstractObjective
    detectionthreshold::DT
    regionlookup::RL
    occurance::OC
    framesperstep::FS
    start::S
end

targets(obj::OffsetRegionObjective) = obj.occurance

simpredictions(obj::OffsetRegionObjective, output) = begin
    regions, steps = size(obj.occurance)
    frames = length(output)
    # Allocate arrays for steps and set all cells to zero 
    outputsteps = [similar(output[1]) for f in 1:steps]
    fill!.(outputsteps, zero(eltype(outputsteps[1])))

    # Get the mean population for steps from the frames in each step:
    # First incomplete step
    firststep = obj.start:obj.framesperstep
    for frame in firststep  
        outputsteps[1] .+= output[frame]
    end
    outputsteps[1] ./= length(firststep) 
    # Full steps
    for frame in obj.framesperstep:frames
        step = stepfromframe(obj.framesperstep, frame)
        outputsteps[step] .+= output[frame]
    end
    # Divide all cells by frames per step to get the mean population
    map(s -> s ./= obj.framesperstep, outputsteps)

    # Allocate a boolean array to contain our presence/absence predictions
    prediction = zeros(Bool, size(obj.occurance))

    # Convert mean cell populations to regional prescence/absence
    for s in 1:steps
        for r in 1:regions
            prediction[r, t] = (sum((obj.regionlookup .== r) .& (outputsteps[s] .> 0)) ./
                                sum((obj.regionlookup .== r))) > obj.detectionthreshold
        end
    end
    prediction 
end


stepfromframe(framesperstep, t) = (t - one(t)) ÷ framesperstep + one(t)

