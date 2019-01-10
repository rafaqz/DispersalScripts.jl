# Precalculate and run the human dispersal model

include("setup.jl")
include("float.jl")

h = h5open("spread_inputs.h5", "r")
names(h)
human = replace(read(h["x_y_popdens"]), NaN=>0)
human = transpose((human[reverse(1:end),:]))
human = human[100:200, 100:200]
RCall.rcall(:image, human.^(1/4))

pg = replace(read(h["x_y_month_intrinsicGrowthRate"]), NaN=>0)
popgrowth = transpose(pg[reverse(1:end),:,1])
popgrowth = popgrowth[100:200, 100:200]
rcall(:image, popgrowth)

dist = distances(popgrowth)
RCall.rcall(:image, (human[1] .* human) .^ (1.7) ./ dist .^ (2e1))

tests = [(human[i, j] == 0 ? 0 : 1) for i = 1:size(human)[1], j = 1:size(human)[2]]
rcall(:image, tests)
# Load Precalc from File
# using JLD2
# @load "human_precalc.jld"

# Run precalc and save to file (uncomment these lines for the first run)
precalc, prop = Dispersal.precalc_human_dispersal(human, 1, 100, 1, 1)
# @save "human_precalc.jld" precalc prop

# Show a single precalc in a Gtk output
init = zeros(size(human))
init[100:110] .= 100.0
single = Dispersal.populate(precalc[50, 50], size(init))
RCall.rcall(:image, Dispersal.populate(precalc[100, 100], size(init)))
# show_frame(GtkOutput(single))

par_a = 1e-5
humandisp = HumanDispersal(precalc, human, par_a)
model = Models(humandisp)
# model = Models(suitability_growth, humandisp)
# model = Models(popdisp, humandisp, suitability_growth)

output = GtkOutput(init, fps=100)
output = ArrayOutput(init, 100)
sim!(output, model, init; tstop=100)
RCall.rcall(:image, output[100])
sum(output[100].>1)
sum(output[100].>1)

resume!(output, model; tadd=1000)
