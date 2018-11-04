# Load packages and raster files

using Revise
using Cellular
using Dispersal
# using CuArrays
# using GPUArrays
# using CUDAnative
using Interact, Blink

using ArchGDAL
function readtiff(file)
    img = ArchGDAL.registerdrivers() do
        ArchGDAL.read(file) do dataset
            ArchGDAL.read(dataset)
        end
    end
    img = img[:,:,1]
end
cropaust(x) = x[3100:3600, 950:1350] # Australia
# cropaust(x) = x[3450:3500, 950:1100] # Queensland

path = "/home/raf/CESAR/Raster/"
growth = cropaust(readtiff(joinpath(path, "new_limited_growth", "limited_growth_2017_01.tif")))
growth_monthly = typeof(growth)[]
for i = 1:9
    push!(growth_monthly, cropaust(readtiff(joinpath(path, "new_limited_growth", "limited_growth_2017_0$i.tif"))))
end
for i = 10:12
    push!(growth_monthly, cropaust(readtiff(joinpath(path, "new_limited_growth", "limited_growth_2017_$i.tif"))))
end
growth_monthly = map(g->exp.(g), growth_monthly)
import Dispersal: pressure 
suit = growth_monthly[1]
suitlay = SuitabilityLayer(suit)
suitseq = SuitabilitySequence((growth_monthly...,), 30);
layers = suitseq

hood = DispersalNeighborhood(; f=exponential, radius=4)
