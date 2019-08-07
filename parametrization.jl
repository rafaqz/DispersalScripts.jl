# Load packages and raster files
# using Distributed
using Pkg: activate
activate(".")
using Revise
using HDF5, Distributed, Statistics, Optim, ColorSchemes, Colors
using CellularAutomataBase, Dispersal, Flatten, LossFunctions, FieldMetadata, Unitful, LabelledArrays
using Unitful: d
using JLD2
Threads.nthreads()
# using DataFrames, CSV # for saving outputs #

# Load simualtions for USA
include("setup_comparison_rulesets.jl")
datafile = "spread_inputs_US_SWD.h5"
my_rulesets, init, tstop, objective = setup_comparison_rulesets(datafile);

# @save "lowlambda.jld2" res
# @save "lowlambda_highres.jld2" res
# @load "lowlambda.jld2" res
# ruleset.rules = reconstruct(ruleset.rules, Optim.minimizer(res))
# ruleset = my_rulesets.full
# @time pars = flatten(ruleset.rules)
# parametriser = Parametriser(ruleset, output, objective, transformfunc, loss, ngroups, groupsize, tstop, threading)
# @time parametriser(pars)
# @profiler parametriser(pars)


# Run the optimizer ######################################################

# @load "SWD/model_comparison_results.jld2" out
# optimresults = @LArray out keys(my_rulesets)

optimresults = @LArray Vector(undef, length(my_rulesets)) keys(my_rulesets)
rulesetkeys = keys(my_rulesets)
rulesetkey = rulesetkeys[2]
transformfunc = x -> 2x - 1
loss = ZeroOneLoss()
threading = Dispersal.DistributedReplicates()
threading = Dispersal.SingleCoreReplicates()
threading = Dispersal.ThreadedReplicates()
groupsize = 5
ngroups = 20
iterations = 1000
output = Dispersal.RegionOutput(init, tstop, objective)
ruleset = my_rulesets[rulesetkey]

for rulesetkey in rulesetkeys
    ruleset = my_rulesets[rulesetkey]
    pars = [flatten(ruleset.rules, Real)...]
    parnames = fieldnameflatten(ruleset.rules, Union{Real,Quantity})
    # Make a labelled array so we can ignore the order
    namedparams = @LArray pars parnames
    show(namedparams)
    parametriser = Parametriser(ruleset, output, objective, transformfunc, loss, ngroups, groupsize, tstop, threading)
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
end

# Accuracy
loss = Dispersal.Accuracy()
accuracy_parametriser = Parametriser(ruleset, output, objective, transformfunc, loss, 
                                     ngroups, groupsize, tstop, threading)
accuracy_parametriser(Optim.minimizer(res))
random = rand(Bool, size(Dispersal.targets(objective)))

# Accuracy of random predictions
LossFunctions.value(loss, Dispersal.(transform.(Dispersal.targets(objective)), transform.(random), AggMode.Sum()))


using DataFrames
paramvectors = LArray{rulesetkeys}([LArray{paramnames}([0.0 for p in paramnames]) for r in rulesetkeys])
sumstats = [:USloss, :USaccuracy, :EUloss, :EUaccuracy]
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
# @load "newresults.jld2" res
include("setup_comparison_rulesets.jl")
datafile = "spread_inputs_US_SWD.h5"
datafile = "spread_inputs_EU_SWD.h5"
datafile = "spread_inputs_Aus_SWD.h5"
datafile = "spread_inputs_Aus_VLM.h5"
my_rulesets, init, tstop, objective = setup_comparison_rulesets(datafile)
ruleset = my_rulesets.full
ruleset.rules = reconstruct(ruleset.rules, Optim.minimizer(res))

for rulesetkey in keys(my_rulesets)
    Parametriser(my_rulesets[rulesetkey], output, objective, transformfunc, loss, ngroups, groupsize, tstop, threading)
    resultdf[rulesetkey][findfirst(x->x==:EUloss, sumstats)] =
        parametriser(Optim.minimizer(optimresults[rulesetkey]))
    parametriser(Optim.minimizer(optimresults[rulesetkey]))
    accuracy_parametriser = Parametriser(ruleset, output, objective, transformfunc, loss, ngroups, groupsize, tstop, threading)
    accuracy_parametriser(Optim.minimizer(res))
end

CellularAutomataBase.showframe(init, output, minimum(init), maximum(init))

CSV.write("SWD/model_comparison.csv", df)


