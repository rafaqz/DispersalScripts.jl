# Run models with floating point arrays - ie for population

include("setup.jl")

minval, maxval = 0.0, 10000.0
init = zeros(Float64, size(init))
init[1:100, 1:500] .= maxval/5

hood = DispersalKernel(; f=exponential, radius=4)
popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=0.1)
pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = [permutedims(pg[:, :, i]) for i in 1:size(pg, 3)]

growth_layers = Sequence(popgrowth, 30);
growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval)
mask_layer = replace(x -> isnan(x) ? 0 : 1, permutedims(read(data["x_y_popdens"])))[:, :, 1]
mask = Dispersal.Mask(mask_layer)
jump = JumpDispersal()

model = Models(popdisp)
model = Models(growth)
model = Models(popdisp, growth)
model = Models(jump, popdisp, growth)
model = Models(jump, (popdisp, growth, mask))
model = Models((popdisp, growth, mask))

# output = ArrayOutput(init, 100)
output = GtkOutput(init, fps=100, store=false, min=minval, max=maxval)
@time sim!(output, model, init; tstop=100000)
maximum(output[100])

resume!(output, model; tadd=10000)
