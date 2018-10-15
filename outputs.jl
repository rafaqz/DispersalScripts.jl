Blink.AtomShell.@dot output.window webContents.setZoomFactor(1.2)
output = Cellular.BlinkOutput(init, model, layers)
# opentools(output.window)

output = Cellular.MuxServer(init, model; port=8000)
output = GtkOutput(init; fps=100, showmax_fps=10) 
output = ArrayOutput(init) 
output = REPLOutput{:braile}(init; fps=800, color=Crayon(foreground=:red, background=:white, bold=true))
output = REPLOutput{:block}(init; fps=50, color=:blue, store=false)
