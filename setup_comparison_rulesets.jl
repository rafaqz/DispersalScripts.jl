# data = h5open("spread_inputs_US.h5", "r")

FloatType = Float32

floatconvert(a::AbstractArray) = convert(Array{FloatType,2}, a)
floatconvert(x::Number) = convert(FloatType, x)

setup_comparison_rulesets(datafile) = begin
    detectionthreshold = floatconvert(1e5)
    month = floatconvert(365.25)d/12
    simtimestep = month
    framesperstep = 12
    data = h5open(datafile, "r")

    # Define parametrization values ########################################
    occurance = convert.(Bool, read(data["state_year_spread"]))
    regionlookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
    steps = size(occurance, 2)
    startmonth = 5 # we start five months/frames in - in May, close to the first sighting.
               # January has strongly negative growth rates in San Jose.
    tstop = steps * framesperstep - startmonth + 1
    objective = OffsetRegionObjective(detectionthreshold, regionlookup, occurance, framesperstep, startmonth)


    # Rules ###########################################################

    # Human
    human_pop = replace(floatconvert.(read(data["x_y_popdens"])), NaN=>missing)
    cellsize = floatconvert(1.0)
    scale = 8
    aggregator = mean
    human_exponent = floatconvert(2.0)
    dist_exponent = floatconvert(2.0)
    dispersalperactivity = floatconvert(5e-9)
    max_dispersers = floatconvert(500.0)
    shortlist_len = 100
    @time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, dispersalperactivity=dispersalperactivity,
                                     max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent,
                                     dist_exponent=dist_exponent, timestep=simtimestep)

    # Constant growth
    constant_growth = ExactLogisticGrowth(intrinsicrate=0.1d^-1)

    # Climate driven growth
    carrycap = floatconvert(100000000.0)
    pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
    popgrowth = [floatconvert(pg[:, :, i]) for i in 1:size(pg, 3)]
    popgrowth = vcat(popgrowth[6:12], popgrowth[1:5])
    # Convert growth arrays to units
    growth_layers = Sequence(popgrowth .* d^-1, month);
    growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=carrycap);

    # Kernel
    λ = floatconvert(1.0)
    radius = 4
    sze = 2radius + 1
    # hood = DispersalKernel{radius}(;kernel= SArray{Tuple{sze,sze}}(zeros(FloatType, sze, sze)) , cellsize=cellsize,
                           # formulation=ExponentialKernel(λ))
    # hood = DispersalKernel{radius}(;kernel = DynamicPaddedArray(zeros(FloatType, sze, sze)),
        # cellsize=cellsize, formulation=ExponentialKernel(λ))

    hood = DispersalKernel{radius}(;kernel=zeros(FloatType, radius, radius), cellsize=cellsize,
                           formulation=ExponentialKernel(λ))
    localdisp = InwardsPopulationDispersal(hood)
    hood.kernel

    # Mask
    masklayer = BitArray(replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"])[:, :, 1]))
    mask = Mask(masklayer)

    # Allee
    minfounders = floatconvert(50.0)
    allee = AlleeExtinction(minfounders=minfounders)

    # Set initialisation array
    init = floatconvert(read(data["x_y_initial"]) .* 10000000) # Arbitrary initial condition
    # masker[:, 300:400] .= false

    # Define combinations for comparison  ##########################
    kwargs = (init=init, timestep=simtimestep, minval=0.0, maxval=carrycap)
    kwargs = (init=init, mask=masklayer, timestep=simtimestep, minval=0.0, maxval=carrycap)
    nolocal = Ruleset(humandisp, (allee, growth); kwargs...)
    noallee = Ruleset(humandisp, (localdisp, growth); kwargs...)
    noclimate = Ruleset(humandisp, (localdisp, allee, constant_growth); kwargs...)
    nohuman = Ruleset((localdisp, allee, growth); kwargs...)
    full = Ruleset(humandisp, (localdisp, allee, growth); kwargs...)
    nolocalnoallee = Ruleset(humandisp, (growth); kwargs...)
    ruleset = Ruleset(humandisp, growth; kwargs...)
    ruleset = Ruleset(growth; kwargs...)
    ruleset = Ruleset(constant_growth; kwargs...)
    ruleset = Ruleset(humandisp; kwargs...)
    ruleset = Ruleset(growth; kwargs...)
    ruleset = Ruleset(humandisp, growth; kwargs...)
    ruleset = Ruleset(localdisp; kwargs...)
    ruleset = Ruleset(localdisp, allee; kwargs...)
    ruleset = Ruleset((localdisp, growth, allee); kwargs...)

    ((full=full, nolocal=nolocal, noallee=noallee, nohuman=nohuman,
      noclimate=noclimate), init, tstop, objective)
      # nolocalnoallee=nolocalnoallee, humangrowth=humangrowth,
      # onlygrowth=onlygrowth, onlyhuman=onlyhuman),
end
