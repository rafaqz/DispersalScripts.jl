using CellularAutomataWeb, Blink

# Frame Processing Colors
truepositivecolor = (0.1, 0.8, 0.2)
falsepositivecolor = (0.8, 0.1, 0.2)
truenegativecolor = (1.0, 1.0, 1.0)
falsenegativecolor = (0.8, 0.8, 0.1)
maskcolor = (0.53, 0.53, 0.53)
regionprocessor = ColorRegionFit(objective, truepositivecolor, falsepositivecolor,
                                truenegativecolor, falsenegativecolor, maskcolor)

# simpleprocessor = GreyscaleZerosProcessor(RGB24(0.5,0.5,0.5))
# schemeprocessor = ColorSchemeProcessor(ColorSchemes.leonardo)
# schemezerosprocessor = ColorSchemeZerosProcessor(ColorSchemes.vermeer, RGB24(0.5, 0.5, 0.5))

output = BlinkOutput(init, ruleset; fps=8, store=true, processor=regionprocessor)#,  extrainit=extrainit,  #summaries=(costs,),
Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0)
