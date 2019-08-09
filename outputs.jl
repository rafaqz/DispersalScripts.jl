output = ArrayOutput(init, tstop)
@time sim!(output, ruleset; tstop=68)

using Crayons
output = REPLOutput(init; style=Block(), fps=10, color=Crayon(foreground=:red, background=:white, bold=true))
output = REPLOutput(init; style=Braile(), fps=5, color=:white, store=false)
output = REPLOutput(init; style=Braile(), fps=10, color=:red, store=false)

truepositivecolor = (0.8, 0.8, 0.8)
falsepositivecolor = (0.8, 0.4, 0.4)
truenegativecolor = (1.0, 1.0, 1.0)
falsenegativecolor = (0.2, 0.3, 0.8)
maskcolor = (0.53, 0.53, 0.53)
processor = ColorRegionFit(objective, truepositivecolor, falsepositivecolor,
                           truenegativecolor, falsenegativecolor, maskcolor)

# Frame Processing Colors
zerocolor = RGB24(0.7) 
maskcolor = RGB24(0.0)
processor = ColorProcessor(Greyscale(), zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.spring, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.Purples_9, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.neon, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.hokusai, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.candy, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.sunset, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.ocean, zerocolor, maskcolor)

processor = ColorProcessor(ColorSchemes.watermelon, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.alpine, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.RdBu_10, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.viridis, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.PRGn_11, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.fall, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.Dark2_3, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.jet, zerocolor, maskcolor)

processor = ColorProcessor(ColorSchemes.Oranges_3, zerocolor, maskcolor)
processor = ColorProcessor(ColorSchemes.PuRd_3, zerocolor, maskcolor)

using CellularAutomataGtk
output = GtkOutput(init .* 0; fps=6, showfps=20, store=true, processor=processor)
output.processor = processor
output.fps = 6 
output.running = false
@time sim!(output, ruleset; tstop=72)

@time sim!(output, ruleset; tstop=72, nreplicates=10)



# using CellularAutomataWeb#, Blink
# output = BlinkOutput(init, ruleset; fps=8, store=true, processor=processor, slider_throttle=1.0)#,  extrainit=extrainit,  #summaries=(costs,),
# Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0)
# output = Cellular.MuxServer(init, model; port=8000)
# output = WebOutput(init, ruleset; fps=8, store=true, processor=processor, slider_throttle=1.0)#,  extrainit=extrainit,  #summaries=(costs,),
# display(output)




@time for i in 1:10 sim!(output, ruleset; tstop=68) end

savegif("sim.gif", output, ruleset; fps=6)


