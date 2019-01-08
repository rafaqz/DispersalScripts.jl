# Load packages and raster files

using Revise
using HDF5 
using Cellular
using Dispersal
using Distributed
using Flatten
 
data = h5open("spread_inputs.h5", "r")
init = read(data["x_y_initial"])
