using Flatten
using Optim

include("setup.jl")

struct Parametriser{M,I,A,Y,S,R,O,CR}
    model::M
    init::I
    args::A
    years::Y
    steps::S
    regions::R
    occurance::O
    region_lookup::CR
end

(p::Parametriser)(a::AbstractVector) = begin
    model = Flatten.reconstruct(p.model, a)
    timesteps = p.years * p.steps_per_year
    s = zeros(Bool, p.regions, p.years)
    output = ArrayOutput(p.init)

    cumsum = 0
    for i = 1:num_runs
        sim!(output, model, p.init, p.args...; time=timesteps)
        for r in 1:p.regions
            for y in 1:p.years
                t = y * p.steps_per_year
                s[r, y] = any(p.region_lookup .== r .&& output[t] .> 0.0)
            end
        end
        cumsum += sum((s .- p.occurance).^2)
    end
    cumsum
end

num_runs = 1000
model = Models(popdisp, humandisp, suitability_growth)
years = 7
regions = 50
steps_per_year = 12

f = Parametriser(model, init, (layers,), years, steps, regions, cell_region, occurance)
optimise(f, flatten(model))
