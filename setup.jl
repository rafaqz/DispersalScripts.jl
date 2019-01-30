# Load packages and raster files

using Revise
using HDF5
using Cellular
using Dispersal
using Distributed
using Flatten
using Unitful: d

data = h5open("spread_inputs_US.h5", "r")
init = read(data["x_y_initial"])
