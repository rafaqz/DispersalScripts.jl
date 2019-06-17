# Load packages and raster files
using Pkg: activate
activate(".")
using Revise, HDF5, Cellular, Dispersal, Distributed, Flatten, Blink, Optim,
      FieldMetadata, Colors, ColorSchemes, LabelledArrays, Statistics, WebIO#, Makie, AbstractPlotting
# @everywhere using Cellular, Dispersal
using Unitful: d
using Dispersal
# data = h5open("spread_inputs_US.h5", "r")
pest_code = "SWD"
data = h5open(string("spread_inputs_US_", pest_code, ".h5"), "r")
names(data)
# Define parametrization values ########################################
minval, maxval = 0.0, 100000.0
month = 365.25d/12
simtimestep = month

# Models ###########################################################
# Human
human_pop = replace(read(data["/x_y_popdens"]), NaN=>missing)
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

# output = SumOutput(init, frames_per_step, steps, Cellular.NullOutput())
# p = RegionParametriser(output, model, init, occurance, region_lookup,
#                  frames_per_step, num_replicates, detection_threshold)
# sim!(output, p.model, p.init; tstop=tstop)

# Get array and field names for the model from Flatten
params = flatten(Vector, model.models)
pnames = fieldnameflatten(model.models)
# Make a labelled array so we can ignore the order
namedparams = @LArray params pnames
# Assign our default parameters to the labelled array
namedparams.minfounders = 23.98
namedparams.param = 0.476303
namedparams.human_exponent = 2.1246
namedparams.dist_exponent = 2.7974
namedparams.par_a = 5.1885e-7
namedparams.max_dispersers = 634.068
show(namedparams)
# Outputs if you want to view/play with the model ######################
model.models = reconstruct(model.models, namedparams)

##################################################################################
# Frame Processing Colors

# Optimisation fit
truepositivecolor = (0.1, 0.1, 0.02)
falsepositivecolor = (0.2, 0.2, 0.2)
truenegativecolor = (1.0, 1.0, 1.0)
falsenegativecolor = (0.8, 0.8, 0.1)
maskcolor = (0.53, 0.53, 0.53)
# regionprocessor = ColorRegionFit(frames_per_step, occurance, region_lookup,
#                            truepositivecolor, falsepositivecolor,
#                            truenegativecolor, falsenegativecolor, maskcolor)
simpleprocessor = GreyscaleZerosProcessor(RGB24(0.5,0.5,0.5))

# Simple colorshceme
schemeprocessor = Cellular.ColorSchemeProcessor(ColorSchemes.leonardo)

# Colorschem with grey zeros
schemezerosprocessor = Cellular.ColorSchemeZerosProcessor(ColorSchemes.vermeer, RGB24(0.5, 0.5, 0.5))

#########################################################################
# Extra Init
# Set these to something useful. The default init will be added to the list.
# extrainit = Dict(:init2 => init, :init3 => init)
#
output = BlinkOutput(init, model; fps=25, store=true, processor=simpleprocessor, # extrainit=extrainit,  #summaries=(costs,),
                      min=minval, max=maxval);
#
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0) # output = REPLOutput{:block}(init)
# save Melbourne incursion
init = read(data["x_y_initial"]) * maxval
output = GtkOutput(init; fps=12, processor=schemeprocessor,min=minval, max=maxval, store=true)
sim!(output, model, init; tstop=12*8)
# savegif("gifs/SWD_melbourne.gif", output)

tstop = 6*12
spreadfileout = string(pest_code, "/spread_sim_save.h5")
rm(spreadfileout)
extrainit = read(data["x_y_initial"])
locnames = collect(keys(extrainit))
reps = 1:10 # simulation replicates
for locname in locnames
      init = extrainit[locname] * maxval
      for rep in reps
            output = ArrayOutput(init, tstop)
            sim!(output, model, init; tstop=tstop)
            out_array = zeros(size(init)...)
            for k = 1:tstop, i= 1:size(init)[1], j = 1:size(init)[2]
                  if out_array[i, j] == 0 && output[k][i, j] > 0
                        out_array[i, j] = k
                  end
            end
            # GtkOutput(out_array[:, :, tstop])
            h5write(spreadfileout, string(locname, "/", rep), out_array)
      end
end

# using Plots
# heatmap(out_array)
#
# dat = h5open("VLM/test2.h5", "r")
# close(dat)
# h5read("VLM/test2.h5", "mygroup2/A")
# h5readattr(spreadfileout, "seisia/1")
# h5open("VLM/output.h5", "w") do file
#     write(file, "seisia", output)  # alternatively, say "@write file A"
# end
