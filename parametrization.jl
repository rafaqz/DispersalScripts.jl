using Flatten
using Optim
using HDF5

include("float.jl")


struct Parametriser{OP,M,I,A,Y,S,R,N,OC,CR}
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
    model.models = Flatten.reconstruct(p.model.models, a)
    timesteps = p.years * p.steps
    s = zeros(Bool, p.regions, p.years)
    cumsum = 0
    for i = 1:p.num_runs
        sim!(p.output, model, p.init, p.args...; time=timesteps)
        for r in 1:p.regions, y in 1:p.years
            t = y * p.steps
            s[r, y] = any((p.region_lookup .== r) .& (output[t] .> 0))
        end
        cumsum += sum((s .- p.occurance).^2)
    end
    cumsum
end

h = h5open("/home/raf/CESAR/spread_inputs.h5", "r")
occurance = convert.(Bool, read(h["state_year_spread"]))
pg = replace(read(h["x_y_month_popgrowthfactor"]), NaN=>0)
popgrowth = [permutedims(pg[:,:,i]) for i in 1:size(pg, 3)]
cell_region = permutedims(convert.(Int, replace(read(h["x_y_state"]), NaN=>0))[:, :, 1])

minmaxrange = 0.0, 10000.0
init = zeros(Float64, size(popgrowth[1]))
init .= minmaxrange[2]
init = ScalableMatrix(init, minmaxrange...)
output = GtkOutput(init, store=true)
model = Models(popdisp, suitability_growth)
# model = Models(popdisp, humandisp, suitability_growth)
years = 6
steps_per_year = 12
timesteps = years * steps_per_year
num_runs = 10
num_regions = maximum(cell_region)
layers = SuitabilitySequence(popgrowth, 1);

p = Parametriser(output, model, init, (layers,), years, steps_per_year, num_regions, num_runs, occurance, cell_region)
p(flatten(model.models))

# o = optimise(p, flatten(Vector, model))
