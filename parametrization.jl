# Load packages and raster files
# using Distributed
using Pkg: activate
activate(".")
using Revise
using HDF5, Distributed, Statistics, Optim, ColorSchemes, Colors
using CellularAutomataBase, Dispersal, Flatten, LossFunctions, FieldMetadata, Unitful, LabelledArrays
using JLD2, DataFrames, CSV
using Unitful: d
Threads.nthreads()
# using DataFrames, CSV # for saving outputs #

# Load simualtions for USA
include("setup_comparison_rulesets.jl")
datafile = "spread_inputs_US_cleaned.h5"
sim_rulesets, init, tstop, objective, output = setup_comparison_rulesets(datafile);
transformfunc = x -> 2x - 1
lossfunc = ZeroOneLoss()

# @save "lowlambda.jld2" res
# @save "lowlambda_highres.jld2" res
rulesetkey = :full
ruleset = sim_rulesets[rulesetkey]
# @load "optimresults2019-08-21T19:27:27.351 .jld2" optimresults
# optimresults
# ruleset.rules = reconstruct(ruleset.rules, Optim.minimizer(optimresults[rulesetkey]))
# @time pars = flatten(ruleset.rules)
# parametriser = Parametriser(ruleset, output, objective, transformfunc, loss, ngroups, groupsize, tstop, threading)
# @time parametriser(pars)
# @profiler parametriser(pars)


# Run the optimizer ######################################################

# @load "SWD/model_comparison_results.jld2" out
# optimresults = @LArray out keys(sim_rulesets)

optimresults = @LArray Vector(undef, length(sim_rulesets)) keys(sim_rulesets)
# threading = Dispersal.DistributedReplicates()
# threading = Dispersal.SingleCoreReplicates()
threading = Dispersal.ThreadedReplicates()
groupsize = 40
ngroups = 5
iterations = 1000

output.running = false
simlosses = @LArray [zeros(groupsize * ngroups) for r in 1:length(sim_rulesets)] keys(sim_rulesets)
for rulesetkey in keys(sim_rulesets)
    println("model: ", rulesetkey)
    ruleset = sim_rulesets[rulesetkey]
    pars = [flatten(ruleset.rules)...]
    ruleset.rules = reconstruct(ruleset.rules, pars .* 0)
    parnames = fieldnameflatten(ruleset.rules, Union{Real,Quantity})
    # Make a labelled array so we can ignore the order
    namedparams = @LArray pars parnames
    show(namedparams)
    parametriser = Parametriser(ruleset, output, objective, transformfunc, lossfunc, ngroups, groupsize, tstop, threading)
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
    simlosses[rulesetkey] = reduce(vcat, parametriser.results)
end


@save "optimresults_$now.jld2" optimresults
@load "optimresults_latest.jld2" optimresults
# optimresults[:full] = res
# Parameter estimates datafram

allparams = Optim.minimizer(optimresults[:full])
paramnames = symbols(allparams)
r = Optim.minimizer(optimresults[:full])
rulesetkeys = keys(sim_rulesets)
paramvectors = LArray{rulesetkeys}([deepcopy(allparams) for r in rulesetkeys]);
paramdf = DataFrame(paramvectors, collect(rulesetkeys))
insert!(paramdf, 1, collect(paramnames), :names)
lims = metaflatten(sim_rulesets[:full].rules, FieldMetadata.limits, Union{Real,Quantity})
lower = [l[1] for l in lims] 
upper = [l[2] for l in lims]
insert!(paramdf, 2, lower, :lowerbound)
insert!(paramdf, 3, upper, :upperbound)

# Loss/accuracy results dataframe
sumstats = [:USloss, :USstd, :USaccuracy, :EUloss, :EUstd, :EUaccuracy]
resultdf = DataFrame(map(k -> k=>zeros(length(sumstats)), keys(optimresults))...)
insert!(resultdf, 1, sumstats, :stat)
accuracy(target, loss) = one(loss) - loss/length(target)

output.running = false
# Loss and accuracy for the USA
for rulesetkey in keys(optimresults)
    println(rulesetkey)
    pars = Optim.minimizer(optimresults[rulesetkey])
    println(pars)
    # Set the parameters to zero, to make sure they are updated in the optimizer
    ruleset = sim_rulesets[rulesetkey]
    ruleset.rules = reconstruct(ruleset.rules, pars .* 0)
    # Make a labelled array so we can ignore the order
    for paramkey in symbols(pars)
        paramkey in paramnames || continue
        paramdf[rulesetkey][paramkey] = pars[paramkey]
    end
    # Fill out the loss for the USA
    parametriser = Parametriser(ruleset, output, objective, transformfunc, 
                                lossfunc, ngroups, groupsize, tstop, threading)
    loss = parametriser(pars)
    resultdf[rulesetkey][findfirst(x->x==:USloss, sumstats)] = loss
    resultdf[rulesetkey][findfirst(x->x==:USstd, sumstats)] = std(vcat(parametriser.results...))
    # Fill out the accuracy for the USA
    resultdf[rulesetkey][findfirst(x->x==:USaccuracy, sumstats)] =
        accuracy(Dispersal.targets(objective), loss)
end


# Loss and accuracy for the EU
include("setup_comparison_rulesets.jl")
datafile = "spread_inputs_EU_cleaned.h5"
# datafile = "spread_inputs_Aus_SWD.h5"
# datafile = "spread_inputs_Aus_VLM.h5"
sim_rulesets, init, tstop, objective, output = setup_comparison_rulesets(datafile)
rulesetkey = :nohuman#:full
ruleset = sim_rulesets[rulesetkey]
ruleset.rules = reconstruct(ruleset.rules, Optim.minimizer(optimresults[rulesetkey]))

for rulesetkey in keys(optimresults)
    println(rulesetkey)
    pars = Optim.minimizer(optimresults[rulesetkey])
    # Set the parameters to zero, to make sure they are updated in the optimizer
    ruleset = sim_rulesets[rulesetkey]
    ruleset.rules = reconstruct(ruleset.rules, pars .* 0) # Fill out the loss for the USA
    parametriser = Parametriser(ruleset, output, objective, transformfunc, 
                                lossfunc, ngroups, groupsize, tstop, threading)
    loss = parametriser(pars)
    resultdf[rulesetkey][findfirst(x->x==:EUloss, sumstats)] = loss
    resultdf[rulesetkey][findfirst(x->x==:EUstd, sumstats)] = std(vcat(parametriser.results...))
    # Fill out the accuracy for the USA
    resultdf[rulesetkey][findfirst(x->x==:EUaccuracy, sumstats)] =
        accuracy(Dispersal.targets(objective), loss)
end

resultdf
paramdf

CSV.write("SWD/model_comparison_results.csv", resultdf)
CSV.write("SWD/model_parameters.csv", paramdf)
