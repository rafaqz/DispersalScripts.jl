# Load packages and raster files

using Revise, HDF5, Cellular, Dispersal, Distributed, Flatten, Blink, Optim, 
      FieldMetadata, Colors, LabelledArrays, Statistics
@everywhere using Cellular, Dispersal
using Unitful: d

data = h5open("spread_inputs_US.h5", "r")
# data = h5open("spread_inputs_EU.h5", "r")
init = read(data["x_y_initial"])



# Define parametrization values ########################################

minval, maxval = 0.0, 100000.0
init .*= maxval
occurance = convert.(Bool, read(data["state_year_spread"]))
region_lookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
steps = size(occurance, 2)
frames_per_step = 12
tstop = steps * frames_per_step
num_replicates = 5
month = 365.25d/12
simtimestep = month
detection_threshold = 0.1



# Models ###########################################################

# Human

human_pop = replace(read(data["x_y_popdens"]), NaN=>missing)
cellsize = 1.0
scale = 8
aggregator = mean
human_exponent = 2.0
dist_exponent = 1.0
par_a = 2.75e-6
max_dispersers = 50.0
shortlist_len = 100
timestep = 1d

@time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, par_a=par_a,
                                 max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent, 
                                 dist_exponent=dist_exponent, timestep=timestep)
# Show a single precalc in a Gtk output
single = Dispersal.populate(humandisp.precalc[20, 20], size(init), scale)
GtkOutput(single, min=0.0, max=maximum(single))
replace(x -> ismissing(x) ? 0.0 : x, humandisp.proportion_covered)

# Growth

pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = [pg[:, :, i] for i in 1:size(pg, 3)]
# Convert growth arrays to units
growth_layers = Sequence(popgrowth * d^-1, month);
growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval);

# Kernel

hood = DispersalKernel(; f=exponential, radius=4)
popdisp = InwardsPopulationDispersal(neighborhood=hood)

# Mask

mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
mask = Dispersal.Mask(mask_layer)

# Allee

allee = AlleeExtinction(minfounders=5.0)



# Define all the possible model combinations  ##########################

model = Models(humandisp; timestep=simtimestep)
model = Models(growth, mask; timestep=simtimestep)
model = Models(humandisp, (growth, mask); timestep=simtimestep)
model = Models((popdisp, growth, mask); timestep=simtimestep)
model = Models(popdisp; timestep=simtimestep)
model = Models((popdisp, allee); timestep=simtimestep)
model = Models(humandisp, (popdisp, growth, mask); timestep=simtimestep)
model = Models((popdisp, allee, growth, mask); timestep=simtimestep)
model = Models(humandisp, (popdisp, allee, growth, mask); timestep=simtimestep);



# Set up the parametrizer ##############################################

output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())
p = RegionParametriser(output, model, init, occurance, region_lookup, 
                 frames_per_step, num_replicates, detection_threshold)
# sim!(output, p.model, p.init; tstop=tstop)

# Get array and field names for the model from Flatten
params = flatten(Vector, model.models)
pnames = fieldnameflatten(model.models)
# Make a labelled array so we can ignore the order
namedparams = LArray{eltype(params),1,pnames}(params)
# Assign our default parameters to the labelled array
namedparams.minfounders = 55.0
namedparams.param = 0.275
namedparams.human_exponent = 1.5
namedparams.dist_exponent = 2.93
namedparams.par_a = 4.5e-7
namedparams.max_dispersers = 112.0
show(namedparams)



# Outputs if you want to view/play with the model ######################
model.models = reconstruct(model.models, namedparams)

truepositivecolor = (0.0, 0.3, 1.0)
falsepositivecolor = (0.5, 0.1, 0.1)
truenegativecolor = (0.0, 0.3, 1.0)
falsenegativecolor = (0.5, 0.1, 0.1)
maskcolor = (0.13, 0.13, 0.13)
processor = ColorRegionFit(frames_per_step, occurance, region_lookup, 
                           truepositivecolor, falsepositivecolor,
                           truenegativecolor, falsenegativecolor, maskcolor)

# output = BlinkOutput(init, model; store=true, processor=processor, min=minval, max=maxval);
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0) # output = REPLOutput{:block}(init)
# sim!(output, model, init; tstop=300)

output = GtkOutput(init; processor=processor, min=minval, max=maxval, store=true)
sim!(output, model, init; tstop=72)

savegif("usa_fit.gif", output)
# @time sim!(output, model, init; tstop=timesteps)
# disp = GtkOutput(init; min=minval, max=maxval*steps, store=false)
# 


# Run the optimizer ######################################################

# Get the lower and upper limits for params with flatten
lims = metaflatten(model.models, FieldMetadata.limits)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]

# p(namedparams)
# o = optimize(p, lower, upper, namedparams, NelderMead())
# res = Optim.optimize(p, lower, upper, namedparams,
                     # SAMIN(), Optim.Options(iterations=100))
