include("setup.jl")

minmaxrange = 0.0, 10000.0
init = ScalableMatrix(zeros(Float64, size(suit)), minmaxrange...)
init[354, 24] = minmaxrange[2]/100

popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.01)
suitability_growth = SuitabilityGrowth(init, 1)

model = Models(suitability_growth)
model = Models(popdisp)
model = Models(popdisp, suitability_growth)

output = GtkOutput(init)
sim!(output, model, init, layers; time=1000)
