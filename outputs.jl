output = ArrayOutput(init, tstop)
@time sim!(output, ruleset; tstop=68)

using Crayons
output = REPLOutput(init; style=Block(), fps=10, color=Crayon(foreground=:red, background=:white, bold=true))
output = REPLOutput(init; style=Braile(), fps=5, color=:white, store=false)
output = REPLOutput(init; style=Braile(), fps=10, color=:red, store=false)

# Frame Processing Colors
processor = GreyscaleProcessor()
processor = GreyscaleZerosProcessor(RGB24(0.5,0.5,0.5))
processor = ColorSchemeProcessor(ColorSchemes.leonardo)
processor = ColorSchemeZerosProcessor(ColorSchemes.vermeer, RGB24(0.5, 0.5, 0.5))

truepositivecolor = (0.1, 0.8, 0.2)
falsepositivecolor = (0.8, 0.1, 0.2)
truenegativecolor = (1.0, 1.0, 1.0)
falsenegativecolor = (0.8, 0.8, 0.1)
maskcolor = (0.53, 0.53, 0.53)
processor = ColorRegionFit(objective, truepositivecolor, falsepositivecolor,
                           truenegativecolor, falsenegativecolor, maskcolor)


using CellularAutomataWeb, Blink
output = BlinkOutput(init, ruleset; fps=8, store=true, processor=processor, slider_throttle=1.0)#,  extrainit=extrainit,  #summaries=(costs,),
Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0)
output = Cellular.MuxServer(init, model; port=8000)
output = WebOutput(init, ruleset; fps=8, store=true, processor=processor, slider_throttle=1.0)#,  extrainit=extrainit,  #summaries=(costs,),
display(output)

using CellularAutomataGtk
output = GtkOutput(init .* 0; fps=6, showfps=20, store=true, processor=processor)
output.fps = 10
output.running = false
@time sim!(output, ruleset; tstop=68)

@time for i in 1:10 sim!(output, ruleset; tstop=68) end

savegif("solocaleu.gif", output, ruleset)
