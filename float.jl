# Run models with floating point arrays - ie for population
include("setup.jl")

popmin = 0.0
popmax = 10000.0
init = zeros(Float64, size(suit))
init[354, 24] = popmax/10
# init .= popmax/10
# init[4, 24] = 10.0 # max/100

popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.5)
carrycap = 100
suitability_growth = SuitabilityLogisticGrowth(suitseq, carrycap)

model = Models(popdisp)
model = Models(suitability_growth)
model = Models(popdisp, suitability_growth)

output = GtkOutput(init, fps=50, store=true, min=popmin, max=popmax)

output = ArrayOutput(init, 100)
@time sim!(output, model, init; tstop=100)
maximum(output[2])
sum(output[2])

resume!(output, model; tadd=100)

output = GtkOutput(suit, fps=50, store=true, min=minimum(suit), max=maximum(suit))

GtkOutput(output[100]/max(output[100]...))
sum(output[])
