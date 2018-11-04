# Choose an output

output = Cellular.BlinkOutput(init, model, layers)
Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.0)
opentools(output.window)

output = Cellular.MuxServer(init, model; port=8000)
output = GtkOutput(init; fps=100, showmax_fps=10) 
output = ArrayOutput(init) 

using Crayons
output = REPLOutput{:braile}(init; fps=800, color=Crayon(foreground=:red, background=:white, bold=true))
output = REPLOutput{:braile}(init; fps=800, color=:red, store=false)
output = REPLOutput{:block}(init; fps=500, color=:blue, store=false)
