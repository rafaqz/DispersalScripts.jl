include("setup.jl")
include("growth.jl")
include("human.jl")

@everywhere using Cellular, Dispersal
using Blink, Optim, FieldMetadata, Colors, LabelledArrays

minval, maxval = 0.0, 100000.0
# init[250:280, 50:80] .= maxval / 100
init .*= maxval
occurance = convert.(Bool, read(data["state_year_spread"]))
region_lookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
frames_per_step = 12
num_replicates = 5
month = 365.25d/12
simtimestep = month

# Convert growth arrays to units
growth_layers = Sequence(popgrowth * d^-1, month);
growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval);

hood = DispersalKernel(; f=exponential, radius=4)
popdisp = InwardsPopulationDispersal(neighborhood=hood)

mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
mask = Dispersal.Mask(mask_layer)

allee = AlleeExtinction(minfounders=5.0)

# Define all the possible models  ######################################
model = Models(humandisp; timestep=simtimestep)
model = Models(growth, mask; timestep=simtimestep)
model = Models(humandisp, (growth, mask); timestep=simtimestep)
model = Models((popdisp, growth, mask); timestep=simtimestep)
model = Models(popdisp; timestep=simtimestep)
model = Models((popdisp, allee); timestep=simtimestep)
model = Models(humandisp, (popdisp, growth, mask); timestep=simtimestep)
model = Models((popdisp, allee, growth, mask); timestep=simtimestep)
model = Models(humandisp, (popdisp, allee, growth, mask); timestep=simtimestep);

# Outputs if you want to view/play with the model ######################

truenegativecolor = maskcolor = RGB24(0.1, 0.1, 0.5)
falsenegativecolor = RGB24(0.5, 0.1, 0.1)
maskcolor = RGB24(0.0, 0.2, 0.2)
processor = ColorRegionFit(frames_per_step, occurance, region_lookup, 
                           truenegativecolor, falsenegativecolor, maskcolor)

# output = BlinkOutput(init, model; store=true, processor=processor, min=minval, max=maxval);
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0) # output = REPLOutput{:block}(init)
# sim!(output, model, init; tstop=300)

# output = GtkOutput(init; processor=processor, min=minval, max=maxval, store=false)
# sim!(output, model, init; tstop=20)

# savegif("usa.gif", output)
# @time sim!(output, model, init; tstop=timesteps)
# disp = GtkOutput(init; min=minval, max=maxval*steps, store=false)

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

# Define parametrization values ########################################

# Assign our default parameters to the labelled array
namedparams.minfounders = 55.0
namedparams.param = 0.275
namedparams.human_exponent = 1.5
namedparams.dist_exponent = 2.93
namedparams.par_a = 4.5e-7
namedparams.max_dispersers = 112.0
show(namedparams)

# Get the lower and upper limits for params with flatten
lims = metaflatten(model.models, FieldMetadata.limits)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]

# Run the optimizer ######################################################
p(namedparams)
# o = optimize(p, lower, upper, namedparams, NelderMead())
# res = Optim.optimize(p, lower, upper, namedparams,
#                      SAMIN(), Optim.Options(iterations=100))
