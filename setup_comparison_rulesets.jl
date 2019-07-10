# data = h5open("spread_inputs_US.h5", "r")

include("objective.jl")

setup_comparison_rulesets(datafile) = begin
    detectionthreshold = 1e5
    month = 365.25d/12
    simtimestep = month
    framesperstep = 12

    data = h5open(datafile, "r")

    # Define parametrization values ########################################
    occurance = convert.(Bool, read(data["state_year_spread"]))
    regionlookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
    steps = size(occurance, 2)
    tstop = steps * framesperstep
    start = 5 # we start five months/frames in - in May, close to the first sighting. 
               # January has strongly negative growth rates in San Jose.
    objective = OffsetRegionObjective(detectionthreshold, regionlookup, occurance, framesperstep, start)


    # Rules ###########################################################

    # Human
    human_pop = replace(read(data["x_y_popdens"]), NaN=>missing)
    cellsize = 1.0
    scale = 8
    aggregator = mean
    human_exponent = 2.0
    dist_exponent = 2.0
    par_a = 5e-8
    max_dispersers = 500.0
    shortlist_len = 100
    @time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, par_a=par_a,
                                     max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent,
                                     dist_exponent=dist_exponent, timestep=simtimestep)

    # Constant growth
    constant_growth = ExactLogisticGrowth(intrinsicrate=0.1d^-1)

    # Climate driven growth
    carrycap = 100000000.0
    pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
    popgrowth = [pg[:, :, i] for i in 1:size(pg, 3)]
    popgrowth = vcat(popgrowth[6:12], popgrowth[1:5])
    # Convert growth arrays to units
    growth_layers = Sequence(popgrowth .* d^-1, month)
    growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=carrycap);

    # Kernel
    param = 0.2
    hood = DispersalKernel(;formulation=ExponentialKernel(param), radius=6)
    localdisp = InwardsPopulationDispersal(neighborhood=hood)

    # Mask
    mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
    mask = Mask(mask_layer)

    # Allee
    minfounders = 100.0 
    allee = AlleeExtinction(minfounders=minfounders)

    # Set initialisation array
    init = read(data["x_y_initial"]) .* 10000000.0 # Arbitrary initial condition 

    # Init Kernel - run one step of spread to account for uncertainy in reporting
    # otherwise we may overfit some parameters to very specific initial conditions
    param = 0.2
    hood = DispersalKernel(;formulation=ExponentialKernel(param), radius=6)
    localdisp = InwardsPopulationDispersal(neighborhood=hood)

    # Define combinations for comparison  ##########################
    kwargs = (init=init, timestep=simtimestep, min=0.0, max=carrycap)
    nolocal = Ruleset(humandisp, (allee, growth, mask); kwargs...)
    noallee = Ruleset(humandisp, (localdisp, growth, mask); kwargs...)
    noclimate = Ruleset(humandisp, (localdisp, allee, constant_growth, mask); kwargs...)
    nohuman = Ruleset((localdisp, allee, growth, mask); kwargs...)
    full = Ruleset(humandisp, (localdisp, allee, growth, mask); kwargs...)
    nolocalnoallee = Ruleset(humandisp, (growth, mask); kwargs...)
    humangrowth = Ruleset(humandisp, growth; kwargs...)
    onlygrowth = Ruleset(growth; kwargs...)
    onlyhuman = Ruleset(humandisp; kwargs...)

    ((full=full, nolocal=nolocal, noallee=noallee, nohuman=nohuman, 
      noclimate=noclimat), init, tstop, objective)
      # nolocalnoallee=nolocalnoallee, humangrowth=humangrowth, 
      # onlygrowth=onlygrowth, onlyhuman=onlyhuman), 
end
