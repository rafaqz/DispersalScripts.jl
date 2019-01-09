# Precalculate and run the human dispersal model

# Load Precalc from File
# using JLD2
# @load "usa_precalc.jld"
human_pop = replace(read(data["x_y_popdens"]), NaN=>missing)
cellsize = 1.0
scale = 8
aggregator = Dispersal.MeanDownsampling()
human_exponent = 4.0/3.0
dist_exponent = 1.0
par_a = 1e-5
shortlist_len = 100


# Run precalc and save to file (uncomment these lines for the first run)
# @time precalc, props = Dispersal.precalc_human_dispersal(downsampled_human, cellsize, shortlist_len, human_exponent, dist_exponent)
# @save "usa_precalc.jld" precalc props

@time humandisp = HumanDispersal(human_pop; scale=scale, shortlist_len=shortlist_len, par_a=par_a,
                                 cellsize=cellsize, human_exponent=human_exponent, dist_exponent=dist_exponent)

# Show a single precalc in a Gtk output
single = Dispersal.populate(humandisp.precalc[20, 20], size(init), scale)
GtkOutput(single, min=0.0, max=maximum(single))
replace(x -> ismissing(x) ? 0.0 : x, humandisp.proportion_covered)
