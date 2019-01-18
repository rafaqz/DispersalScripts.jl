include("setup.jl")
include("growth.jl")
include("human.jl")
# using Blink
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

@everywhere (p::Parametriser)(params) = begin
    # Rebuild the model with the current parameters
    names = fieldnameflatten(p.model.models) 
    println("Parameters: ", collect(zip(names, params)))
    p.model.models = Flatten.reconstruct(p.model.models, params)
    # Calculate the timespan
    tstop = p.years * p.steps
    s = zeros(Bool, p.regions, p.years)
    # Used a labelled array so we don't have to know the index number
    # fort par_a.
    lv = LVector{eltype(params),typeof(params),names}(params)
    # Large par a take forever. Limit the max value.
    cumsum = @distributed (+) for i = 1:p.num_runs
        o = deepcopy(p.output)
        sim!(o, p.model, p.init; tstop=tstop)
        for t in 1:p.years
            for r in 1:p.regions
                s[r, t] = any((p.region_lookup .== r) .& (o[t] .> 0))
            end
        end
        val = sum((s .- p.occurance).^2)
        println("replicate: ", i, " - result: ", val)
        val
    end
    println("output: ", cumsum, "\n")
    cumsum
end

minval, maxval = 0.0, 100000.0
init[250:280, 50:80] .= maxval / 100
occurance = convert.(Bool, read(data["state_year_spread"]))
cell_region = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
years = 6
steps_per_year = 12
tstop = years * steps_per_year
num_runs = 15
num_regions = maximum(cell_region)
month = 365.25d/12
simtimestep = month

hood = DispersalKernel(; f=exponential, radius=4, param=param)
popdisp = InwardsPopulationDispersal(neighborhood=hood)

growth_layers = Sequence(popgrowth * d^-1, month); 
growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval);

mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
mask = Dispersal.Mask(mask_layer)

allee = AlleeExtinction(minfounders=5)

model = Models(humandisp; timestep=simtimestep)
model = Models(growth, mask; timestep=simtimestep)
model = Models(humandisp, (growth, mask); timestep=simtimestep)
model = Models((popdisp, growth, mask); timestep=simtimestep)
model = Models(popdisp; timestep=simtimestep)
model = Models((popdisp, allee); timestep=simtimestep)
model = Models(humandisp, (popdisp, growth, mask); timestep=simtimestep)
model = Models((popdisp, allee, growth, mask); timestep=simtimestep)
model = Models(humandisp, (popdisp, allee, growth, mask); timestep=simtimestep);

# output = BlinkOutput(init, model; min=minval, max=maxval);
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0) # output = REPLOutput{:block}(init)
# sim!(output, model, init; tstop=300)

# output = GtkOutput(init; min=minval, max=maxval, store=false)
# sim!(output, model, init; tstop=200)

# savegif("usa.gif", output)
# @time sim!(output, model, init; tstop=timesteps)
# disp = GtkOutput(init; min=minval, max=maxval*steps_per_year, store=false)

output = SumOutput(init, steps_per_year, years, Cellular.NullOutput())
p = Parametriser(output, model, init, years, steps_per_year, num_regions,
                 num_runs, occurance, cell_region)
# sim!(output, p.model, p.init; tstop=tstop)

minfounders = 140.0
param = 1.0

params = flatten(Vector, model.models)
reconstruct(model.models, params) 
names = fieldnameflatten(Vector, model.models)

lims = metaflatten(model.models, FieldMetadata.limits)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]

# @time p(params)
o = optimize(p, lower, upper, params)
