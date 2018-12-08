# Precalculate and run the human dispersal model

include("setup.jl")
include("float.jl")

world = readtiff(joinpath(path, "population_density.tif"))
human = max.(0.0, cropaust(world)) # .^(4/3)
humanlay = HumanLayer(human)

# Load Precalc from File
using JLD2
@load "human_precalc.jld"

# Run precalc and save to file (uncomment these lines for the first run)
# precalc, prop = Dispersal.precalc_human_dispersal(human, 1, 100)
# @save "human_precalc.jld" precalc prop

# Show a single precalc in a Gtk output
single = Dispersal.populate(precalc[50, 50], size(init))
show_frame(GtkOutput(single))

humandisp = HumanDispersal(precalc=precalc, prob_threshold=0.5)
model = Models(humandisp)
model = Models(suitability_growth, humandisp)
model = Models(popdisp, humandisp, suitability_growth)

output = GtkOutput(init, fps=100)
sim!(output, model, init; tstop=1000)

resume!(output, model; tadd=1000)
