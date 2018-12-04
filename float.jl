# Run models with floating point arrays - ie for population

include("setup.jl")

minmaxrange = 0.0, 10000.0
init = ScalableMatrix(zeros(Float64, size(suit)), minmaxrange...)
init[354, 24] = minmaxrange[2]/100
# init[4, 24] = 10.0 # minmaxrange[2]/100

popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.0001)
suitability_growth = SuitabilityExponentialGrowth(init)

model = Models(suitability_growth)
model = Models(popdisp)
# model = Models(popdisp, suitability_growth)

# output = GtkOutput(init, fps=300, store=true)

output = ArrayOutput(init)
@time sim!(output, model, init, layers; tstop=500)

resume!(output, model, layers; tadd=100)
