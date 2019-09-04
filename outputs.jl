output = ArrayOutput(init, tstop)
@time sim!(output, ruleset; tstop=68)

using Crayons
output = REPLOutput(init; style=Block(), fps=10, color=Crayon(foreground=:red, background=:white, bold=true))
output = REPLOutput(init; style=Braile(), fps=5, color=:white, store=false)
output = REPLOutput(init; style=Braile(), fps=10, color=:red, store=false)


scheme = Greyscale()
scheme = ColorSchemes.spring
scheme = ColorSchemes.Purples_9
scheme = ColorSchemes.neon
scheme = ColorSchemes.hokusai
scheme = ColorSchemes.candy
scheme = ColorSchemes.sunset
scheme = ColorSchemes.ocean

scheme = ColorSchemes.watermelon
scheme = ColorSchemes.alpine
scheme = ColorSchemes.RdBu_10
scheme = ColorSchemes.viridis
scheme = ColorSchemes.PRGn_11
scheme = ColorSchemes.fall
scheme = ColorSchemes.Dark2_3
scheme = ColorSchemes.jet

scheme = ColorSchemes.Oranges_3
scheme = ColorSchemes.PuRd_3

# Frame Processing Colors
zerocolor = RGB24(0.7) 
maskcolor = RGB24(0.0)
processor = ColorProcessor(scheme, zerocolor, maskcolor)

truescheme = Greyscale(max=0.4, min=0.0)
falsescheme = ColorSchemes.RdGy_5
falsescheme = ColorSchemes.autumn1
falsezerocolor = RGB24(0.80, .85, .27)
truezerocolor = RGB24(1.0, 1.0, 1.0)
maskcolor = RGB24(0.53, 0.53, 0.53)
processor = ColorRegionFit(objective, truescheme, falsescheme, 
                           truezerocolor, falsezerocolor, maskcolor)

rulesetkey = :full
rulesetkey = :noclimate
rulesetkey = :nolocal
rulesetkey = :noallee
rulesetkey = :nohuman
ruleset = sim_rulesets[rulesetkey]
ruleset.rules = reconstruct(ruleset.rules, Optim.minimizer(optimresults[rulesetkey]))

using CellularAutomataGtk
output = GtkOutput(init .* 0; fps=6, showfps=20, store=true, processor=processor)
output.fps = 10
output.running = false
output.processor = processor
@time sim!(output, ruleset; tstop=103)

@time sim!(output, ruleset; tstop=68, nreplicates=10)


using CellularAutomataWeb#, Blink
output = BlinkOutput(init, ruleset; fps=8, store=true, processor=processor, slider_throttle=1.0)#,  extrainit=extrainit,  #summaries=(costs,),
output.interface.running = false
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0)
# output = Cellular.MuxServer(init, model; port=8000)
# output = WebOutput(init, ruleset; fps=8, store=true, processor=processor, slider_throttle=1.0)#,  extrainit=extrainit,  #summaries=(costs,),
# display(output)


savegif("sim.gif", output, ruleset; fps=6)


