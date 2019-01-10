using Revise
using Cellular
using Dispersal
using Interact
using RCall
using HDF5

grid_dims = (100, 100)
tmax = 1000

function plotgrowth(output)
    # function for plotting pop size vs time in R
    x = [i for i =1:tmax]
    y = [sum(output[i]) for i = 1:tmax]
    RCall.rcall(:plot, x, y, "l", ylab = "y")
end

### simple local dispersal simulation that works###
init = convert(Array{Float64}, zeros(grid_dims))
init[(Int.(round.(grid_dims)./2))...] = 100000.0
hood = DispersalKernel(; f=exponential, radius=2)
localdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=1.0)
model = Models(localdisp)
output = ArrayOutput(init, tmax)
@time sim!(output, model, init; tstop=tmax)
plotgrowth(output)
RCall.rcall(:image, output[700])

# make exponential population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
exp_growth = EulerExponentialGrowth(intrinsicrate = 0.01, timestep = 1)
model = Models(exp_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
plotgrowth(output)

# make logistic population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
logistic_growth = EulerLogisticGrowth(intrinsicrate = 0.1, timestep = 1, carrycap = 100)
model = Models(logistic_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
plotgrowth(output)

# make exact logistic population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
ex_log_growth = ExactLogisticGrowth(intrinsicrate = 0.1, timestep = 1, carrycap = 100)
model = Models(ex_log_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
plotgrowth(output)

# suitability growth layer and logistic population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
popgrowth = fill(0.1, (grid_dims..., 10))
popgrowth = [permutedims(popgrowth[:,:,i]) for i in 1:size(popgrowth, 3)]
for i = 6:10
    popgrowth[i] .= -0.1
end
popgrowthseq = Sequence((popgrowth...,), 5);
suit_log_growth = SuitabilityEulerLogisticGrowth(layers = popgrowthseq, carrycap = 100)
model = Models(suit_log_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
plotgrowth(output)
RCall.rcall(:image, output[end])

### local dispersal + suitabilty logistic growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[(Int.(round.(grid_dims)./2))...] = 1
model = Models(localdisp, ex_log_growth)
output = ArrayOutput(init, tmax)
@time sim!(output, model, init; tstop = tmax)
plotgrowth(output)
RCall.rcall(:image, output[50])

### simulation with real seasonalith data for US ###
# load data
h = h5open("spread_inputs.h5", "r")
names(h)
pg = replace(read(h["x_y_month_intrinsicGrowthRate"]), NaN=>0)
human = replace(read(h["x_y_popdens"]), NaN=>0)
@time RCall.rcall(:image, transpose((human[reverse(1:end),:]).^(1/4)))
popgrowth = [permutedims(pg[reverse(1:end),:,i]) for i in 1:size(pg, 3)]
RCall.rcall(:image, popgrowth[5])

init = zeros(Float64, size(popgrowth[1]))
init[304, 24] = 10000
tmax = 365
suitseq = Sequence((popgrowth...,), 30);
suit_ex_log_growth = SuitabilityExactLogisticGrowth(layers = suitseq, carrycap = 100)
suit_log_growth = SuitabilityEulerLogisticGrowth(layers = suitseq, carrycap = 100)
hood = DispersalKernel(; f=exponential, radius=4)
popdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=1)
allee = AlleeExtinction(minfounders = 1)
model = Models(localdisp, suit_log_growth, allee)
output = ArrayOutput(init, tmax)
@time sim!(output, model, init; tstop=tmax)
plotgrowth(output)
RCall.rcall(:image, output[end])

## Gtk doesnt work
# output = GtkOutput(init, fps=50, store=false, min=0, max=100)
# @time sim!(output, model, init; tstop=tmax)

# resume!(output, model, layers; time=4000
