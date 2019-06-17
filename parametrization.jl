# Load packages and raster files

using Pkg: activate
activate(".")
using Revise, HDF5, Mux, Cellular, Dispersal, Distributed, Flatten, InteractBase, InteractBulma, Blink, Optim, LossFunctions,
      FieldMetadata, Colors, ColorSchemes, LabelledArrays, Statistics, WebIO#, Makie, AbstractPlotting
# @everywhere using Cellular, Dispersal
using Unitful: d
using DataFrames, Optim, LabelledArrays, CSV, JLD2 # for saving outputs

# data = h5open("spread_inputs_US.h5", "r")
include("setup_comparison_models.jl")
my_models = setup_comparison_models("spread_inputs_US_SWD.h5");
model = my_models.model_full;


################## REMOVE DEPENDENCY OF THIS STUFF FOR BLINK AND DELETE #######
data = h5open("spread_inputs_US_SWD.h5", "r")
minval, maxval = 0.0, 100000.0
init = read(data["x_y_initial"])
init .*= maxval
occurance = convert.(Bool, read(data["state_year_spread"]))
region_lookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
steps = size(occurance, 2)
frames_per_step = 12
tstop = steps * frames_per_step
loss = ZeroOneLoss()
month = 365.25d/12
simtimestep = month
num_replicates = 1
detection_threshold = 0.1

##################################################################################
# Frame Processing Colors

# Optimisation fit
truepositivecolor = (0.1, 0.1, 0.02)
falsepositivecolor = (0.2, 0.2, 0.2)
truenegativecolor = (1.0, 1.0, 1.0)
falsenegativecolor = (0.8, 0.8, 0.1)
maskcolor = (0.53, 0.53, 0.53)
regionprocessor = ColorRegionFit(frames_per_step, occurance, region_lookup,
                           truepositivecolor, falsepositivecolor,
                           truenegativecolor, falsenegativecolor, maskcolor)
simpleprocessor = GreyscaleZerosProcessor(RGB24(0.5,0.5,0.5))

# Simple colorshceme
schemeprocessor = Cellular.ColorSchemeProcessor(ColorSchemes.leonardo)

# Colorschem with grey zeros
schemezerosprocessor = Cellular.ColorSchemeZerosProcessor(ColorSchemes.vermeer, RGB24(0.5, 0.5, 0.5))

########################### BLINK OUTPUT #############################
# output = BlinkOutput(init, model; fps=25, store=true, processor=regionprocessor, # extrainit=extrainit,  #summaries=(costs,),
#                       min=minval, max=maxval)
#
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0) # output =

# Run the optimizer ######################################################

# Get the lower and upper limits for params with flatten
lims = metaflatten(model.models, FieldMetadata.limits)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]

############## compare models ########################
params = flatten(Vector, model.models)
pnames = fieldnameflatten(model.models)
namedparams = @LArray params pnames

objective = RegionObjective(detection_threshold, region_lookup, occurance)

output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())
p = Parametriser(output, model, init, num_replicates,
                 objective, loss, tstop)
p(namedparams)

# my_models = (model_noclimate,)
@load "SWD/model_comparison_results.jld2" out
# model = my_models[3]
for model in my_models
      # Get array and field names for the model from Flatten
      params = flatten(Vector, model.models)
      pnames = fieldnameflatten(model.models)
      # Make a labelled array so we can ignore the order
      namedparams = @LArray params pnames
      show(namedparams)
      # Outputs if you want to view/play with the model ######################
      model.models = reconstruct(model.models, namedparams)

      output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())
      p = Parametriser(output, model, init, num_replicates,
                       objective, loss, tstop)

      # Get the lower and upper limits for params with flatten
      lims = metaflatten(model.models, FieldMetadata.limits)
      lower = [l[1] for l in lims]
      upper = [l[2] for l in lims]
      #
      # p(namedparams)
      res = Optim.optimize(p, lower, upper, namedparams,
                          SAMIN(), Optim.Options(iterations=3))
      push!(out, res)
end
@save "SWD/model_comparison_results.jld2" out

# make dictionary to match up parameters from simpler models
syms = push!(collect(symbols(Optim.minimizer(out[1]))), :USloss,:USaccuracy, :EUloss, :EUaccuracy)
fvars = Dict()
for i in 1:length(syms)
   fvars[syms[i]] = i
end

# fill table with paramters and loss function values
df = DataFrame(row = syms; map(m -> (m = 0.0), my_models)...)
for i in 1:length(my_models)
      ind = [fvars[key] for key in collect(symbols(Optim.minimizer(out[i])))]
      df[i + 1][ind] = convert(Array, Optim.minimizer(out[i]))
      df[i + 1][fvars[:USloss]] = Optim.minimum(out[i])
end
df

# fill loss for EU
my_models = setup_comparison_models("spread_inputs_EU_SWD.h5");
model = my_models.model_full;
data = h5open("spread_inputs_EU_SWD.h5", "r")
minval, maxval = 0.0, 100000.0
init = read(data["x_y_initial"])
init .*= maxval
objective = convert.(Bool, read(data["state_year_spread"]))
region_lookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
steps = size(objective, 2)
frames_per_step = 12
tstop = steps * frames_per_step
output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())
for i in 1:length(my_models)
      p = Parametriser(output, my_models[i], init, objective, region_lookup, frames_per_step, num_replicates, detection_threshold)
      df[i + 1][fvars[:EUloss]] = p(Optim.minimizer(out[i]))
end
CSV.write("SWD/model_comparison.csv", df)

Optim.minimizer(out[4])
