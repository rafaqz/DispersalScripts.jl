# data = h5open("spread_inputs_US.h5", "r")
setup_comparison_models(datafile) = begin
    data = h5open(datafile, "r")
    init = read(data["x_y_initial"])

    # Define parametrization values ########################################
    minval, maxval = 0.0, 100000.0
    init .*= maxval
    occurance = convert.(Bool, read(data["state_year_spread"]))
    region_lookup = convert.(Int, replace(read(data["x_y_state"]), NaN=>0))[:, :, 1]
    steps = size(occurance, 2)
    frames_per_step = 12
    tstop = steps * frames_per_step
    month = 365.25d/12
    simtimestep = month
    detection_threshold = 0.1

    # Models ###########################################################

    # Human
    minfounders = 23.98

    human_pop = replace(read(data["x_y_popdens"]), NaN=>missing)
    cellsize = 1.0
    scale = 8
    aggregator = mean
    param = 0.476303
    human_exponent = 2.1246
    dist_exponent = 2.7974
    par_a = 5.1885e-7
    max_dispersers = 634.068
    minfounders = 23.98
    shortlist_len = 100
    timestep = 1d # this is bug!!!
    @time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, par_a=par_a,
                                     max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent,
                                     dist_exponent=dist_exponent, timestep=timestep)


    # Constant growth
    constant_growth = ExactLogisticGrowth(intrinsicrate=0.1d^-1)
    # Climate driven growth
    pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
    popgrowth = [pg[:, :, i] for i in 1:size(pg, 3)]

    # Convert growth arrays to units
    growth_layers = Sequence(popgrowth * d^-1, month);
    growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=maxval);

    # Kernel
    hood = DispersalKernel(;param=param, f=exponential, radius=4)
    popdisp = InwardsPopulationDispersal(neighborhood=hood)
    # Mask
    mask_layer = replace(x -> isnan(x) ? 0 : 1, read(data["x_y_popdens"]))[:, :, 1]
    mask = Dispersal.Mask(mask_layer)

    # Allee
    allee = AlleeExtinction(minfounders=minfounders)

    # Define combinations for comparison  ##########################
    model_spread = Models(humandisp, (allee, growth, mask); timestep=simtimestep);
    model_noallee = Models(humandisp, (popdisp, growth, mask); timestep=simtimestep);
    model_noclimate = Models(humandisp, (popdisp, allee, constant_growth, mask); timestep=simtimestep);
    model_nohuman = Models((popdisp, allee, growth, mask); timestep=simtimestep);
    model_full = Models(humandisp, (popdisp, allee, growth, mask); timestep=simtimestep);

    (model_full=model_full, model_noallee=model_noallee, model_nohuman=model_nohuman, model_noclimate=model_noclimate)

end
