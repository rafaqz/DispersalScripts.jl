# Load packages and raster files

using Revise
using Cellular
using Dispersal
using Interact, Blink, HDF5 
 
h = h5open("spread_inputs.h5", "r")
occurance = convert.(Bool, read(h["state_year_spread"]))
pg = replace(read(h["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = [permutedims(pg[:,:,i]) for i in 1:size(pg, 3)]
cell_region = permutedims(convert.(Int, replace(read(h["x_y_state"]), NaN=>-2))[:, :, 1])

suit = popgrowth[1]
suitlay = popgrowth[1]
suitseq = Sequence((popgrowth...,), 30);

hood = DispersalKernel(; f=exponential, radius=4)


