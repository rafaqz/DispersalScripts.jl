include("setup.jl")

using Flatten, Optim, HDF5, Distributed


@everywhere using Cellular

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
        sim!(p.output, model, p.init, p.args...; tstop=timesteps)
        for r in 1:p.regions, y in 1:p.years
            t = y * p.steps
            annual_presence = reduce(+, a)
            s[r, y] = any((p.region_lookup .== r) .& (annual_presence .> 0))
        end
        out = sum((s .- p.occurance).^2)
        println(out)
        out
    end
    cumsum
end

h = h5open("spread_inputs.h5", "r")
occurance = convert.(Bool, read(h["state_year_spread"]))
pg = replace(read(h["x_y_month_popgrowthfactor"]), NaN=>0)
popgrowth = [permutedims(pg[:,:,i]) for i in 1:size(pg, 3)]
cell_region = permutedims(convert.(Int, replace(read(h["x_y_state"]), NaN=>0))[:, :, 1])

minmaxrange = 0.0, 10000.0
init = zeros(Float64, size(popgrowth[1]))
init .= minmaxrange[2]
init = ScalableMatrix(init, minmaxrange...)
# output = GtkOutput(init; fps=5000, store=true)
output = ArrayOutput(init)

popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.0001)
suitability_growth = SuitabilityExponentialGrowth(init)

model = Models(popdisp, suitability_growth)
layers = SuitabilitySequence(popgrowth, 1);
# model = Models(popdisp, humandisp, suitability_growth)
years = 6
steps_per_year = 12
timesteps = years * steps_per_year
num_runs = 1
num_regions = maximum(cell_region)

p = Parametriser(output, model, init, (layers,), years, steps_per_year, num_regions, num_runs, occurance, cell_region)
using Profile, ProfileView
Profile.clear()
@profile p(flatten(model.models))
ProfileView.view()

# o = optimize(p, flatten(Vector, model))
