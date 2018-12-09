# Run models with floating point arrays - ie for population
# james comment
include("setup.jl")

min = 0.0 
max = 10000.0
init = zeros(Float64, size(suit))
init[354, 24] = max/10
init .= max/10
# init[4, 24] = 10.0 # max/100

popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.1)
suitability_growth = SuitabilityExponentialGrowth(suitseq, min, max)

model = Models(popdisp)
model = Models(suitability_growth)
model = Models(popdisp, suitability_growth)

output = GtkOutput(init, fps=50, store=true, min=min, max=max)

output = ArrayOutput(init, 100)
@time sim!(output, model, init; tstop=100)
maximum(output[100])

resume!(output, model; tadd=100)

output = GtkOutput(suit, fps=50, store=true, min=minimum(suit), max=maximum(suit))

