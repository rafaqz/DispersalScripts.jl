# Load packages and raster files

using Pkg: activate
activate(".")
using Revise, HDF5, CellularAutomataBase, Dispersal, Distributed, Flatten, LossFunctions,
      FieldMetadata, Colors, Unitful, ColorSchemes, LabelledArrays, Statistics, Optim, LabelledArrays
using Unitful: d
# using DataFrames, CSV, JLD2 # for saving outputs # 

# Constants for all simulations
loss = ZeroOneLoss()
nreplicates = 1

include("setup_comparison_rulesets.jl")
my_rulesets, init, tstop, objective = setup_comparison_rulesets("spread_inputs_US_SWD.h5");

##################################################################################
# Frame Processing Colors

# # Optimisation fit
# truepositivecolor = (0.1, 0.1, 0.02)
# falsepositivecolor = (0.2, 0.2, 0.2)
# truenegativecolor = (1.0, 1.0, 1.0)
# falsenegativecolor = (0.8, 0.8, 0.1)
# maskcolor = (0.53, 0.53, 0.53)
# regionprocessor = ColorRegionFit(objective, truepositivecolor, falsepositivecolor,
#                                  truenegativecolor, falsenegativecolor, maskcolor)
# simpleprocessor = GreyscaleZerosProcessor(RGB24(0.5,0.5,0.5))

# # Simple colorshceme
# schemeprocessor = ColorSchemeProcessor(ColorSchemes.leonardo)

# # Colorschem with grey zeros
# schemezerosprocessor = ColorSchemeZerosProcessor(ColorSchemes.vermeer, RGB24(0.5, 0.5, 0.5))

########################### BLINK OUTPUT #############################
# output = BlinkOutput(init, ruleset; fps=25, store=true, processor=regionprocessor)#,  extrainit=extrainit,  #summaries=(costs,),
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0) # output =

# Run the optimizer ######################################################

# @load "SWD/model_comparison_results.jld2" out
# optimresults = @LArray out keys(my_rulesets)

optimresults = @LArray Vector(undef, length(my_rulesets)) keys(my_rulesets)

sumstats = [:USloss,:USaccuracy, :EUloss, :EUaccuracy]
paramnames = union(fieldnameflatten.(getfield.(values(my_rulesets), :rules))...)
rulesetnames = keys(my_rulesets)
rulesetname = rulesetnames[1]
f = my_rulesets.No_climate.rules[2][3]
ruleset = reconstruct(f, flatten(f, Real), Real)


# ruleset = my_ruleset[3]
for rulesetname in rulesetnames
    ruleset = my_rulesets[rulesetname] 
    # Get array and field names for the ruleset from Flatten
    params = [flatten(ruleset.rules, Real)...]
    # Use Quantity for fieldnames to get unitful field names instead of :val
    pnames = fieldnameflatten(Foo(1u"m", 2.0, 3u"g"), Union{Real,Quantity})
    # Make a labelled array so we can ignore the order
    namedparams = @LArray params pnames
    show(namedparams)
    parametriser = Parametriser(ruleset, objective, loss, nreplicates, tstop)
    # Get the lower and upper limits for params with flatten
    lims = metaflatten(ruleset.rules, FieldMetadata.limits, Union{Real,Quantity})
    lower = [l[1] for l in lims]
    upper = [l[2] for l in lims]
    res = Optim.optimize(parametriser, lower, upper, namedparams, SAMIN(), Optim.Options(iterations=3))
    optimresults[rulesetname] = res
end


fieldnameflatten(my_rulesets[:No_climate].rules[2][3])

paramvectors = LArray{rulesetnames}([LArray{paramnames}([0.0 for p in paramnames]) for r in rulesetnames])
resultdf = DataFrame(row = sumstats; map(m -> (m = 0.0), my_rulesets)...)

for rulesetkey in keys(my_rulesets)
    params = Optim.minimizer(optimresults[rulesetkey])
    lims = metaflatten(my_rulesets[rulesetkey].rules, FieldMetadata.limits)
    for paramkey in symbols(params)
        paramvectors[rulesetkey][paramkey] = params[paramkey]
    end
    resultdf[rulesetkey][findfirst(x->x==:USloss, sumstats)] = 
        Optim.minimum(optimresults[rulesetkey])
end

paramdf = DataFrame(paramvectors, collect(rulesetnames))
insert!(paramdf, 1, collect(paramnames), :names)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]
insert!(paramdf, 2, lower, :lower)
insert!(paramdf, 3, upper, :upper)
paramdf


# fill loss for EU
my_rulesets, init, tstop, objective = setup_comparison_rulesets("spread_inputs_EU_SWD.h5");
for rulesetkey in keys(my_rulesets)
    parametriser = Parametriser(my_rulesets[rulesetkey], objective, loss, nreplicates, tstop)
    resultdf[rulesetkey][findfirst(x->x==:EUloss, sumstats)] = 
        parametriser(Optim.minimizer(optimresults[rulesetkey]))
end

CSV.write("SWD/model_comparison.csv", df)

Optim.minimizer(optimresults[4])
