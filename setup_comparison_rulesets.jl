# data = h5open("spread_inputs_US.h5", "r")

import FieldMetadata: @relimits, limits
import Flatten: @reflattenable flatten


# Human
@relimits struct HumanDispersal
    human_exponent  | (0.0, 3.0)
    dist_exponent   | (0.0, 3.0)
    dispersalperpop | (0.0, 1e-8)
    max_dispersers  | (1e1, 1e5)
end
# Constant growth
@relimits struct ExactLogisticGrowth
    intrinsicrate | (0.0, 10.0)
end
# Dont parametrise carrycap
@reflattenable struct SuitabilityExactLogisticGrowth
    carrycap | false
end
# Kernel
@relimits struct ExponentialKernel
    λ | (0.0, 0.8)
end
# Alee
@relimits struct AlleeExtinction
    minfounders | (1e0, 1e3)
end

FloatType = Float32

floatconvert(a::AbstractArray) = convert(Array{FloatType,2}, a)
floatconvert(x::Number) = convert(FloatType, x)

setup_comparison_rulesets(datafile) = begin
    data = h5open(datafile, "r")
    cellsize = floatconvert(1.0)
    init = Dispersal.downsample(floatconvert(read(data["x_y_initial"]) .* 10000000), mean, Int(cellsize)) # Arbitrary initial condition
    masklayer = BitArray(replace(x -> isnan(x) ? 0 : 1, Dispersal.downsample(read(data["x_y_popdens"])[:, :, 1], mean, Int(cellsize))))
    month = floatconvert(365.25)d/12
    simtimestep = month
    framesperstep = 12

    # Define parametrization objective ########################################
    detectionthreshold = floatconvert(1e7)
    occurance = convert.(Bool, read(data["state_year_spread"]))
    regionlookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
    steps = size(occurance, 2)
    startmonth = 5 # we start five months/frames in - in May, close to the first sighting.
               # January has strongly negative growth rates in San Jose.
    tstop = steps * framesperstep - startmonth + 1
    objective = RegionObjective(detectionthreshold, regionlookup, occurance, framesperstep, startmonth)


    # Rules ###########################################################


    human_pop = replace(floatconvert.(read(data["x_y_popdens"])), NaN=>missing)
    display(human_pop[100:300,100:300])
    scale = 8
    aggregator = mean
    human_exponent = floatconvert(2.0)
    dist_exponent = floatconvert(2.0)
    dispersalperpop = floatconvert(1e-9)
    max_dispersers = floatconvert(500.0)
    shortlist_len = 100
    @time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, dispersalperpop=dispersalperpop,
                                     max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent,
                                     dist_exponent=dist_exponent, timestep=simtimestep)

    constant_growth = ExactLogisticGrowth(intrinsicrate=0.1d^-1)

    # Climate driven growth
    carrycap = floatconvert(1e8)
    pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
    popgrowth = [Dispersal.downsample(floatconvert(pg[:, :, i]), mean, Int(cellsize)) for i in 1:size(pg, 3)]
    popgrowth = vcat(popgrowth[6:12], popgrowth[1:5])
    # Convert growth arrays to units
    growth_layers = Sequence(popgrowth .* d^-1, month);
    growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=carrycap);

    λ = floatconvert(0.05)
    radius = 3
    sze = 2radius + 1
    hood = DispersalKernel{radius}(;kernel=zeros(FloatType, radius, radius), cellsize=cellsize,
                           formulation=ExponentialKernel(λ))
    localdisp = InwardsPopulationDispersal(hood)
    display(hood.kernel .* carrycap)

    minfounders = floatconvert(5.0)
    allee = AlleeExtinction(minfounders=minfounders)


    # Define combinations for comparison  ##########################
    kwargs = (init=init, mask=masklayer, timestep=simtimestep, minval=0.0, maxval=carrycap)
    nolocal = Ruleset(humandisp, (allee, growth); kwargs...)
    noallee = Ruleset(humandisp, (localdisp, growth); kwargs...)
    noclimate = Ruleset(humandisp, (localdisp, allee, constant_growth); kwargs...)
    nohuman = Ruleset((localdisp, allee, growth); kwargs...)
    full = Ruleset(humandisp, (localdisp, allee, growth); kwargs...)

    # Manual
    kwargs = (init=init, timestep=simtimestep, minval=0.0, maxval=carrycap)
    ruleset = Ruleset(growth; kwargs...)
    ruleset = Ruleset(constant_growth; kwargs...)
    ruleset = Ruleset(humandisp; kwargs...)
    ruleset = Ruleset(growth; kwargs...)
    ruleset = Ruleset(humandisp, growth; kwargs...)
    ruleset = Ruleset(localdisp; kwargs...)
    ruleset = Ruleset(localdisp, allee; kwargs...)
    ruleset = Ruleset(humandisp, (growth, allee); kwargs...)
    ruleset = Ruleset((localdisp, growth, allee); kwargs...)

    ((full=full, nolocal=nolocal, noallee=noallee, nohuman=nohuman,
      noclimate=noclimate), init, tstop, objective)
      # nolocalnoallee=nolocalnoallee, humangrowth=humangrowth,
      # onlygrowth=onlygrowth, onlyhuman=onlyhuman),
end
