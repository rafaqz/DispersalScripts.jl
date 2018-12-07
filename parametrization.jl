include("setup.jl")

@everywhere using Cellular
@everywhere import Cellular: store_frame!

using Flatten, Optim, Distributed

@Ok @Frames @everywhere struct SumOutput{SF} <: AbstractArrayOutput{T} 
    sum_frames::SF
end

@everywhere SumOutput(frames::AbstractVector, sum_frames, num_frames) = begin
    o = SumOutput{typeof(frames), typeof(sum_frames)}(frames, [false], sum_frames)
    allocate!(o, frames[1], 2:num_frames)
    o
end

@everywhere store_frame!(o::SumOutput, frame, t) = begin
    sze = size(o[1])
    ts = (t + o.sum_frames - 2one(t)) ÷ o.sum_frames + 1
    for j in 1:sze[2]
        for i in 1:sze[1]
            @inbounds o[ts][i, j] += frame[i, j]
        end
    end
end

@everywhere struct Parametriser{OP,M,I,A,Y,S,R,N,OC,CR}
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

@everywhere (p::Parametriser)(a) = begin
    model.models = Flatten.reconstruct(p.model.models, a)
    tstop = p.years * p.steps
    s = zeros(Bool, p.regions, p.years)
    cumsum = @distributed (+) for i = 1:p.num_runs
        o = deepcopy(p.output)
        sim!(o, model, p.init, p.args...; tstop=timesteps)
        for t in 1:p.years
            for r in 1:p.regions 
                s[r, t] = any((p.region_lookup .== r) .& (o[t] .> 0))
            end
        end
        out = sum((s .- p.occurance).^2)
        println(out)
        out
    end
    cumsum
end

h = h5open("spread_inputs.h5", "r")
occurance = convert.(Bool, read(h["state_year_spread"]))
pg = replace(read(h["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = [permutedims(pg[:,:,i]) for i in 1:size(pg, 3)]
cell_region = permutedims(convert.(Int, replace(read(h["x_y_state"]), NaN=>0))[:, :, 1])

minmaxrange = 0.0, 10000.0
init = zeros(Float64, size(popgrowth[1]))
init[200:300,200:300] .= minmaxrange[2]
init = ScalableMatrix(init, minmaxrange...)
# output = GtkOutput(init; fps=5000, store=true)

popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.0001)
suitability_growth = SuitabilityExponentialGrowth(init)

model = Models(popdisp, suitability_growth)
layers = SuitabilitySequence(popgrowth, 1);
# model = Models(popdisp, humandisp, suitability_growth)
years = 6
steps_per_year = 12
timesteps = years * steps_per_year + 1
num_runs = 100
num_regions = maximum(cell_region)
output = SumOutput(init, steps_per_year, years + 1)

p = Parametriser(output, model, init, (layers,), years, steps_per_year, num_regions, num_runs, occurance, cell_region)
@time p(flatten(model.models))

# using Profile, ProfileView
# Profile.clear()
# @profile p(flatten(model.models))
# ProfileView.view()

# o = optimize(p, flatten(Vector, model))
