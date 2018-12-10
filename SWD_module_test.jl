using Revise
using Cellular
using Dispersal
using Interact
using RCall

grid_dims = (100, 100)
tmax = 1000

### simple local dispersal simulation that works###
init = convert(Array{Float64}, zeros(grid_dims))
init[50,50] = 100.0
hood = DispersalKernel(; f=exponential, radius=2)
localdisp = InwardsPopulationDispersal(neighborhood=hood, fraction=1.0)
model = Models(localdisp)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
# test if the total number of individuals is constant through time
[sum(output[i]) for i = 1:10]
RCall.rcall(:image, output[50])


# make exponential population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
exp_growth = ExponentialGrowth(growthrate = 0.01)
model = Models(exp_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
x = [i for i =1:tmax]
y = [sum(output[i]) for i = 1:tmax]
RCall.rcall(:plot, x, y, "l", ylab = "y")


# make logistic population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
logistic_growth = LogisticGrowth(growthrate = 0.1, carrycap = 100)
model = Models(logistic_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
x = [i for i =1:tmax]
y = [sum(output[i]) for i = 1:tmax]
RCall.rcall(:plot, x, y, "l", ylab = "y")



# suitability growth layer and logistic population growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[1,1] = 10
popgrowth = fill(0.1, (grid_dims..., 10))
popgrowth = [permutedims(popgrowth[:,:,i]) for i in 1:size(popgrowth, 3)]
# for i = 6:10
#     popgrowth[i] .= -0.1
# end
popgrowthseq = Sequence((popgrowth...,), 5);
carrycap = 100
suit_log_growth = SuitabilityLogisticGrowth(popgrowthseq, carrycap)
model = Models(suit_log_growth)
tmax = 100
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
output
x = [i for i =1:tmax]
y = [sum(output[i]) for i = 1:tmax]
RCall.rcall(:plot, x, y, "l", ylab = "y")
RCall.rcall(:image, output[1])

### local dispersal + suitabilty logistic growth simulation
init = convert(Array{Float64}, zeros(grid_dims))
init[(Int.(round.(grid_dims)./2))...] = 10
model = Models(localdisp, suit_log_growth)
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
x = [i for i =1:tmax]
y = [sum(output[i]) for i = 1:tmax]
RCall.rcall(:plot, x, y, "l", ylab = "y")
RCall.rcall(:image, output[5])


### simple simulation that works###
h = h5open("spread_inputs.h5", "r")
pg = replace(read(h["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = [permutedims(pg[reverse(1:end),:,i]) for i in 1:size(pg, 3)]
RCall.rcall(:image, popgrowth[5])
suitseq = Sequence((popgrowth...,), 30);
suit_log_growth = SuitabilityLogisticGrowth(suitseq, 1000)
model = Models(popdisp, suitability_growth)
init = zeros(Float64, size(popgrowth[1]))
init[304, 24] = popmax

tmax = 100
output = ArrayOutput(init, tmax)
sim!(output, model, init; tstop=tmax)
x = [i for i = 1:tmax]
y = [sum(output[i]) for i = 1:tmax]
RCall.rcall(:plot, x, y, "l", ylab = "y")
RCall.rcall(:image, output[5])

pos1 = popgrowth[1]
pos2 = [pos1[i, j] > 0 ? 1 : 0 for i in 1:679, j in 1:509]
RCall.rcall(:image,  pos2)

# resume!(output, model, layers; time=4000)
