include("setup.jl")
include("growth.jl")
include("human.jl")
# using Blink
# using Optim 

@everywhere using Cellular
@everywhere import Cellular: store_frame!, show_frame, allocate_frames!, @Ok, @Frames 


@Ok @Frames @everywhere struct SumOutput{SF} <: AbstractArrayOutput{T} 
    steps::SF
end

@everywhere show_frame(output::SumOutput, t::Number) = nothing

# An output that sums frames for a number of frames
@everywhere SumOutput(frames::AbstractVector, steps, years) = begin
    o = SumOutput{typeof(frames), typeof(steps)}(frames, [false], steps)
    allocate_frames!(o, frames[1], 2:years)
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
end

@everywhere year_from_month(o, t) = (t - 1one(t)) ÷ o.steps + one(t) @everywhere struct Parametriser{OP,M,I,A,Y,S,R,N,OC,CR}
    output::OP
    model::M
    init::I
    args::A
    years::Y
    steps::S
    regions::R
    num_runs::N
    occurance::OC
    region_lookup::CR
end

(p::Parametriser)(a) = begin
    # Rebuild the model with the current parameters
    p.model.models = Flatten.reconstruct(p.model.models, a)
    # Calculate the timespan
    tstop = p.years * p.steps
    s = zeros(Bool, p.regions, p.years)
    cumsum = @distributed (+) for i = 1:p.num_runs
        o = deepcopy(p.output)
        sim!(o, p.model, p.init, p.args...; tstop=tstop)
        for t in 1:p.years
            for r in 1:p.regions 
                s[r, t] = any((p.region_lookup .== r) .& (o[t] .> 0))
            end
        end
        out = sum((s .- p.occurance).^2)
        out
    end
    cumsum
end

minval, maxval = 0.0, 100000.0
init[350:470, 200:320] .= maxval / 100
occurance = convert.(Bool, read(data["state_year_spread"]))
cell_region = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
years = 6
steps_per_year = 12
tstop = years * steps_per_year
num_runs = 12
num_regions = maximum(cell_region)

hood = DispersalKernel(; f=exponential, radius=4)
popdisp = InwardsPopulationDispersal(neighborhood=hood)
growth_layers = Sequence(popgrowth, 30);
growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval);
mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
allee = AlleeExtinction(100)
mask = Dispersal.Mask(mask_layer)
model = Models(humandisp)
model = Models(growth, mask)
model = Models(humandisp, (growth, mask))
model = Models((popdisp, growth, mask))
model = Models(humandisp, (popdisp, growth, mask))
model = Models((popdisp, growth, allee, mask))
model = Models(humandisp, (popdisp, growth, allee, mask))

output = BlinkOutput(init, model; min=minval, max=maxval)
Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0)
#
# output = REPLOutput{:block}(init)
# sim!(output, model, init; tstop=300)

# output = GtkOutput(init; min=minval, max=maxval, store=false)
# sim!(output, model, init; tstop=200)
# savegif("usa.gif", output)
# @time sim!(output, model, init; tstop=timesteps)

# output = SumOutput(init, steps_per_year, years)

# p = Parametriser(output, model, init, (layers,), years, steps_per_year, num_regions, 
                 # num_runs, occurance, cell_region)
# @time p(flatten(model.models))

# output = ArrayOutput(init, 100)
# using Profile, ProfileView
# Profile.clear()
# resume!(output, model; tadd=500)
# ProfileView.view()

# @profile p(flatten(model.models))
# o = optimize(p, flatten(Vector, model))
