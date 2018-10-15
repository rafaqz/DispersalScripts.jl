include("setup.jl")

world = readtiff(joinpath(path, "population_density.tif"))
human = cropaust(world)
humanlay = HumanLayer(human)

# Run precalc
human .^= 4/3
precalc, props = Dispersal.precalc_human_dispersal(human, 1, 200)
precalc[100,140]

# Show a single precalc in a Gtk output
out = GtkOutput(human)
fill!(out[1], 0.0)
Dispersal.populate!(out[1], precalc[300,150])
show_frame(out)

humandisp = HumanDispersal(precalc=precalc)

model = Models(humandisp)
output = GtkOutput(init)

sim!(output, model, init, layers; time=100)
