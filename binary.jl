# Run models on binary 1/0 arrays

include("setup.jl")

init = zeros(Int64, size(growth))
init[354, 24] = 1
# init = CuArray(init)

localdisp = InwardsBinaryDispersal(neighborhood=hood)
jumpdisp = JumpDispersal()
suitability_mask = SuitabilityMask()

model = Models(jumpdisp)
model = Models(localdisp)
model = Models(suitability_mask)
model = Models(jumpdisp, localdisp)
model = Models(localdisp, suitability_mask)
model = Models(jumpdisp, suitability_mask)
model = Models(localdisp, suitability_mask, jumpdisp)
model = Models(localdisp, jumpdisp, suitability_mask)

output = GtkOutput(init, fps=100)
sim!(output, model, init, layers; time=4000)
resume!(output, model, layers; time=4000)
