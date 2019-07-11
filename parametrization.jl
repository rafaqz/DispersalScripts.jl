# Load packages and raster files

@everywhere using Pkg: activate
@everywhere activate(".")
@everywhere using HDF5, Distributed, Statistics, Optim, ColorSchemes, Colors 
@everywhere using CellularAutomataBase, Dispersal, Flatten, LossFunctions,
      FieldMetadata, Unitful, LabelledArrays
@everywhere using Unitful: d
# using DataFrames, CSV, JLD2 # for saving outputs # 

# Constants for all simulations

@everywhere include("setup_comparison_rulesets.jl")
my_rulesets, init, tstop, objective = setup_comparison_rulesets("spread_inputs_US_SWD.h5");
ruleset = my_rulesets.full
ruleset = my_rulesets.noallee
ruleset = my_rulesets.nohuman
ruleset = my_rulesets.nolocal
ruleset = my_rulesets.noclimate

# Run the optimizer ######################################################

# @load "SWD/model_comparison_results.jld2" out
# optimresults = @LArray out keys(my_rulesets)

optimresults = @LArray Vector(undef, length(my_rulesets)) keys(my_rulesets)
sumstats = [:USloss,:USaccuracy, :EUloss, :EUaccuracy]
paramnames = (union(fieldnameflatten.(getfield.(values(my_rulesets), :rules))...)...,)
rulesetnames = keys(my_rulesets)
rulesetname = rulesetnames[4]
transform = x -> 2x - 1
loss = ZeroOneLoss()
nreplicates = 10
iterations = 10


for rulesetname in rulesetnames
    ruleset = my_rulesets[rulesetname] 
    # Get array and field names for the ruleset from Flatten
    params = [flatten(ruleset.rules, Real)...]
    # Use Quantity for fieldnames to get unitful field names instead of :val
    pnames = fieldnameflatten(ruleset.rules, Union{Real,Quantity})
    # Make a labelled array so we can ignore the order
    namedparams = @LArray params pnames
    show(namedparams)
    parametriser = Parametriser(ruleset, objective, transform, loss, nreplicates, tstop)
    # Get the lower and upper limits for params with flatten
    lims = metaflatten(ruleset.rules, FieldMetadata.limits, Union{Real,Quantity})
    lower = [l[1] for l in lims]
    upper = [l[2] for l in lims]
    @time res = Optim.optimize(parametriser, lower, upper, namedparams, SAMIN(), 
                               Optim.Options(iterations=iterations,
                                             show_trace=true,
                                             store_trace=true
                                            ))
    optimresults[rulesetname] = res
end



@save "optimresults.jld2" optimresults

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
