# Load packages and raster files
# using Distributed
using Pkg: activate
activate(".")
using Revise
using HDF5, Distributed, Statistics, Optim, ColorSchemes, Colors
using CellularAutomataBase, Dispersal, Flatten, LossFunctions, FieldMetadata, Unitful, LabelledArrays
using Unitful: d
using JLD2, FileIO
Threads.nthreads()
# using DataFrames, CSV # for saving outputs #

# Constants for all simulations
include("objective.jl")
include("setup_comparison_rulesets.jl")
datafile = "spread_inputs_US_SWD.h5"
my_rulesets, init, tstop, objective = setup_comparison_rulesets(datafile);
ruleset = my_rulesets.full

ruleset.rules[1].precalc
ruleset = my_rulesets.nohuman
ruleset = my_rulesets.noallee
ruleset = my_rulesets.nolocal
ruleset = my_rulesets.noclimate

# @save "fullmodel_results.jld2" res
@load "fullmodel_results.jld2" res
@load "full2.jld2" res
ruleset.rules = reconstruct(ruleset.rules, Optim.minimizer(res))

# Run the optimizer ######################################################

# @load "SWD/model_comparison_results.jld2" out
# optimresults = @LArray out keys(my_rulesets)

optimresults = @LArray Vector(undef, length(my_rulesets)) keys(my_rulesets)
sumstats = [:USloss, :USaccuracy, :EUloss, :EUaccuracy]
paramnames = (union(fieldnameflatten.(getfield.(values(my_rulesets), :rules))...)...,)
rulesetkeys = keys(my_rulesets)
rulesetkey = rulesetkeys[1]
transformfunc = x -> 2x - 1
loss = ZeroOneLoss()
nreplicates = 500
iterations = 5

# for rulesetname in rulesetkeys
ruleset = my_rulesets[rulesetkey]
# Get array and field names for the ruleset from Flatten
pars = [flatten(ruleset.rules, Real)...]
# Use Quantity for fieldnames to get unitful field names instead of :val
pnames = fieldnameflatten(ruleset.rules, Union{Real,Quantity})
# Make a labelled array so we can ignore the order
namedparams = @LArray pars pnames
show(namedparams)
parametriser = Parametriser(ruleset, objective, transformfunc, loss, nreplicates, tstop)
# Get the lower and upper limits for params with flatten
lims = metaflatten(ruleset.rules, FieldMetadata.limits, Union{Real,Quantity})
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]
res = Optim.optimize(parametriser, lower, upper, namedparams, SAMIN(),
                     Optim.Options(iterations=iterations,
                                   show_trace=true,
                                   store_trace=true
                                  ))

optimresults[rulesetkey] = res


loss = Dispersal.Accuracy()
parametriser = Parametriser(ruleset, objective, transformfunc, loss, nreplicates, tstop)
parametriser(flatten(ruleset.rules))



using DataFrames
paramvectors = LArray{rulesetkeys}([LArray{paramnames}([0.0 for p in paramnames]) for r in rulesetkeys])
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

paramdf = DataFrame(paramvectors, collect(rulesetkeys))
insert!(paramdf, 1, collect(paramnames), :names)
lower = [l[1] for l in lims]
upper = [l[2] for l in lims]
insert!(paramdf, 2, lower, :lower)
insert!(paramdf, 3, upper, :upper)
paramdf


# fill loss for EU
datafile = "spread_inputs_EU_SWD.h5"
my_rulesets, init, tstop, objective = setup_comparison_rulesets(datafile)
for rulesetkey in keys(my_rulesets)
    parametriser = Parametriser(my_rulesets[rulesetkey], objective, transformfunc, loss, nreplicates, tstop)
    resultdf[rulesetkey][findfirst(x->x==:EUloss, sumstats)] =
        parametriser(Optim.minimizer(optimresults[rulesetkey]))
    parametriser(Optim.minimizer(optimresults[rulesetkey]))
end

CSV.write("SWD/model_comparison.csv", df)
