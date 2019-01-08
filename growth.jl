pg = replace(read(data["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = [pg[:, :, i] for i in 1:size(pg, 3)]

growth_layers = Sequence(popgrowth, 30)
suitability_growth = SuitabilityExactLogisticGrowth(layers=growth_layers, carrycap=10)
