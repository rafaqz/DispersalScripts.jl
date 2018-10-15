using Revise
using Cellular
using Dispersal
# using Crayons
# using CuArrays
# using GPUArrays
# using CUDAnative
# using JLD2
using Interact, Blink # using Blink # using Mux # using Flatten
# using CLArrays
# using FileIO 
# using ImageView 

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

# cropaust(x) = x[950:1250, 3400:3600] # Australia
path = "/home/raf/CESAR/Raster/"
growth = readtiff(joinpath(path, "new_limited_growth", "limited_growth_2017_01.tif"))
growth_monthly = typeof(growth)[]
for i = 1:9
    push!(growth_monthly, cropaust(readtiff(joinpath(path, "new_limited_growth", "limited_growth_2017_0$i.tif"))))
end
for i = 10:12
    push!(growth_monthly, cropaust(readtiff(joinpath(path, "new_limited_growth", "limited_growth_2017_$i.tif"))))
end
import Dispersal: pressure 
suit = growth_monthly[1] 
suitlay = SuitabilityLayer(growth_monthly[1])
suitseq = SuitabilitySequence((growth_monthly...,), 30);
layers = suitseq

init = zeros(Int64, size(growth))
init[354, 24] = 1
# init = CuArray(init)
hood = DispersalNeighborhood(; f=exponential, radius=4, init=init)
