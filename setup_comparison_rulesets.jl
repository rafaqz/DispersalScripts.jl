# data = h5open("spread_inputs_US.h5", "r")
#
setup_comparison_rulesets(datafile) = begin
    detectionthreshold = 0.1
    month = 365.25d/12
    simtimestep = month
    framesperstep = 12

    data = h5open(datafile, "r")

    # Define parametrization values ########################################
    occurance = convert.(Bool, read(data["state_year_spread"]))
    regionlookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
    steps = size(occurance, 2)
    tstop = steps * framesperstep
    objective = RegionObjective(detectionthreshold, regionlookup, occurance, framesperstep)


    # Rules ###########################################################

    # Human
    human_pop = replace(read(data["x_y_popdens"]), NaN=>missing)
    cellsize = 1.0
    scale = 8
    aggregator = mean
    human_exponent = 2.1246
    dist_exponent = 2.7974
    par_a = 5.1885e-7
    max_dispersers = 634.068
    shortlist_len = 100
    @time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, par_a=par_a,
                                     max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent,
                                     dist_exponent=dist_exponent, timestep=simtimestep)

    # Constant growth
    constant_growth = ExactLogisticGrowth(intrinsicrate=0.1d^-1)

    # Climate driven growth
    carrycap = 1000000.0
    pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
    popgrowth = [pg[:, :, i] for i in 1:size(pg, 3)]
    # Convert growth arrays to units
    growth_layers = Sequence(popgrowth * d^-1, simtimestep)
    growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=carrycap);
    flattenable(ExactLogisticGrowth())

    # Kernel
    param = 0.476303
    hood = DispersalKernel(;formulation=ExponentialKernel(param), radius=4)
    localdisp = InwardsPopulationDispersal(neighborhood=hood)

    # Mask
    mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
    mask = Mask(mask_layer)

    # Allee
    minfounders = 23.98
    allee = AlleeExtinction(minfounders=minfounders)

    # Set initialisation array
    init = read(data["x_y_initial"]) .* 10000.0 # Arbitrary initial condition 

    # Define combinations for comparison  ##########################
    kwargs = (init=init, timestep=simtimestep, min=0.0, max=carrycap)
    ruleset_nolocal = Ruleset(humandisp, (allee, growth, mask); kwargs...)
    ruleset_noallee = Ruleset(humandisp, (localdisp, growth, mask); kwargs...)
    ruleset_noclimate = Ruleset(humandisp, (localdisp, allee, constant_growth, mask); kwargs...)
    ruleset_nohuman = Ruleset((localdisp, allee, growth, mask); kwargs...)
    ruleset_full = Ruleset(humandisp, (localdisp, allee, growth, mask); kwargs...)

    ((Full=ruleset_full, No_local=ruleset_nolocal, No_allee=ruleset_noallee, 
     No_human=ruleset_nohuman, No_climate=ruleset_noclimate), 
     init, tstop, objective)
end
