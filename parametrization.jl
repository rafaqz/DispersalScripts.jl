using Flatten
using Optim

struct Parametriser{M,Y,S,R,O}
    model::M
    years::Y
    steps::S
    regions::R
    occurance::O
end

(p::Parametriser)(a) = begin
    model = Flatten.reconstruct(p.model, a)
    timesteps = p.years * p.steps_per_year
    s = zeros(Bool, p.regions, p.years)
    output = ArrayOutput(init)

    for i = 1:num_runs
        sim!(output, model, init, layers; time = timesteps)
        for r in 1:p.regions
            for y in 1:p.years
                s[r, y] = p.occurance .== r .&& output[y * p.steps_per_year] .> 0.0
            end
        end
    end
    sum((s .- p.occurance).^2)
end

include("setup.jl")

num_runs = 1000
model = ModelList(popdisp, humandisp, suitability_growth)
years = 7
regions = 50
steps_per_year = 12

f = Parametriser(model, years, steps, regions, occurance)
optimise(f, flatten(model))
