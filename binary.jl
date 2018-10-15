# Run models on binary 1/0 arrays

include("setup.jl")

init = zeros(Int64, size(growth))
init[354, 24] = 1
# init = CuArray(init)

localdisp = InwardsBinaryDispersal(neighborhood=hood)
suitability_mask = SuitabilityMask()
jumpdisp = JumpDispersal()
model = Models(localdisp, suitability_mask)
model = Models(jumpdisp, localdisp)
model = Models(localdisp)
model = Models(jumpdisp)

output = GtkOutput(init)
sim!(output, model, init, layers; time=200)
