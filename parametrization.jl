include("setup.jl")
using Blink
using Optim
using FieldMetadata

@everywhere using LabelledArrays
@everywhere using Cellular
@everywhere import Cellular: store_frame!, show_frame, allocate_frames!, @Ok, @Frames, is_async
@everywhere import Flatten: flatten


@Ok @Frames @everywhere struct SumOutput{DI,SF} <: AbstractArrayOutput{T} where DI <: AbstractOutput disp::DI
    steps::SF
end

@everywhere show_frame(output::SumOutput, t::Number) = nothing

# An output that sums frames for a number of frames
@everywhere SumOutput(frames::AbstractVector, steps::Number, years::Number, disp::AbstractOutput) = begin
    o = SumOutput{typeof.((frames, disp, steps))...}(frames, [false], disp, steps)
    allocate_frames!(o, frames[1], 2:years)
    map(f -> f .= 0, o.frames)
    o
end
# Sums frames on the fly to reduce storage
@everywhere store_frame!(o::SumOutput, frame, t) = begin
    sze = size(o[1])
    # Determine the timestep being summed to
    ts = year_from_month(o, t)
    # Add frame to current sum frame
    for j in 1:sze[2]
        for i in 1:sze[1]
            @inbounds o[ts][i, j] += frame[i, j]
        end
    end
    show_frame(o.disp, o[ts], t)
end

@everywhere year_from_month(o, t) = (t - 1one(t)) ÷ o.steps + one(t)

@everywhere struct Parametriser{OP,M,I,Y,S,R,N,OC,CR}
    output::OP
    model::M
    init::I
    years::Y
    steps::S
    regions::R
    num_runs::N
    occurance::OC
    region_lookup::CR
end

@everywhere crossentropy = function(y, p, minprob = 1e-9)
    p = min.(p, 1 - minprob)
    p = max.(p, minprob)
    -sum( y  .* log.(p) .+ (ones(size(y)) .- y) .* log.(ones(size(p)) .- p))
end

@everywhere (p::Parametriser)(params) = begin
    # Rebuild the model with the current parameters
    names = fieldnameflatten(p.model.models)
    println("Parameters: ", collect(zip(names, params)))
    p.model.models = Flatten.reconstruct(p.model.models, params)
    # Calculate the timespan
    tstop = p.years * p.steps
    s = zeros(Bool, p.regions, p.years)
    det_thresh = 0.1 # proportion of state infected for detection
    # Used a labelled array so we don't have to know the index number
    # fort par_a.
    # Large par a take forever. Limit the max value.
    # function to pass back s array (region x year)
    cumsum = @distributed (+) for i = 1:p.num_runs
        o = deepcopy(p.output)
        sim!(o, p.model, p.init; tstop=tstop)
        for t in 1:p.years
            for r in 1:p.regions

                s[r, t] = (Base.sum((p.region_lookup .== r) .& (o[t] .> 0)) ./
                                Base.sum((p.region_lookup .== 6) )) > det_thresh
            end
        end
        val = sum((s .== p.occurance)) / prod(size(p.occurance))
        println("replicate: ", i, " - accuracy: ", val)
        s
    end
    # then build and array of s array means
    probs = cumsum ./ p.num_runs
    loss = crossentropy(p.occurance, probs)
    println("cross-entropy loss: ", loss, "\n")
    loss
end

# Define parametrization values ########################################

minval, maxval = 0.0, 100000.0
# init[250:280, 50:80] .= maxval / 100
init .*= maxval
occurance = convert.(Bool, read(data["state_year_spread"]))
cell_region = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
years = size(data["state_year_spread"]))[2]
steps_per_year = 12
tstop = years * steps_per_year
num_runs = 5
num_regions = maximum(cell_region)
month = 365.25d/12
simtimestep = month

# Define model components ##############################################
include("human.jl")

include("growth.jl")
# Convert growth arrays to units
growth_layers = Sequence(popgrowth * d^-1, month);
growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval);

hood = DispersalKernel(; f=exponential, radius=4)
popdisp = InwardsPopulationDispersal(neighborhood=hood)

mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
mask = Dispersal.Mask(mask_layer)

allee = AlleeExtinction(minfounders=5.0)

# Define all the possible models  ######################################
model = Models(humandisp; timestep=simtimestep)
model = Models(growth, mask; timestep=simtimestep)
model = Models(humandisp, (growth, mask); timestep=simtimestep)
model = Models((popdisp, growth, mask); timestep=simtimestep)
model = Models(popdisp; timestep=simtimestep)
model = Models((popdisp, allee); timestep=simtimestep)
model = Models(humandisp, (popdisp, growth, mask); timestep=simtimestep)
model = Models((popdisp, allee, growth, mask); timestep=simtimestep)
model = Models(humandisp, (popdisp, allee, growth, mask); timestep=simtimestep);

# Outputs if you want to view/play with the model ######################

# output = BlinkOutput(init, model; min=minval, max=maxval);
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.5) # output = REPLOutput{:block}(init)
# sim!(output, model, init; tstop=300)

# output = GtkOutput(init; min=minval, max=maxval, store=false)
# sim!(output, model, init; tstop=24)

# savegif("usa.gif", output)
# @time sim!(output, model, init; tstop=timesteps)
# disp = GtkOutput(init; min=minval, max=maxval*steps_per_year, store=false)

# Set up the parametrizer ##############################################

output = SumOutput(init, steps_per_year, years, Cellular.NullOutput())
p = Parametriser(output, model, init, years, steps_per_year, num_regions,
                 num_runs, occurance, cell_region)
# sim!(output, p.model, p.init; tstop=tstop)

# Get array and field names for the model from Flatten
params = flatten(Vector, model.models)
pnames = fieldnameflatten(model.models)
# Make a labelled array so we can ignore the order
namedparams = LArray{eltype(params),1,pnames}(params)

# Assign our default parameters to the labelled array
namedparams.minfounders = 55.0
namedparams.param = 0.275
namedparams.human_exponent = 1.5
namedparams.dist_exponent = 2.93
namedparams.par_a = 4.5e-7
namedparams.max_dispersers = 112.0
show(namedparams)

# Get the lower and upper limits for params with flatten
lims = metaflatten(model.models, FieldMetadata.limits)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]

# Run the optimizer ######################################################
p(namedparams)
# o = optimize(p, lower, upper, namedparams, NelderMead())
# res = Optim.optimize(p, lower, upper, namedparams,
#                      SAMIN(), Optim.Options(iterations=100))
