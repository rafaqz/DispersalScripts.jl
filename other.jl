
output = GtkOutput(init, fps=100, store=true)
sim!(output, model, init, layers; time=1000)
savegif("output.gif", output)

using PerceptualColourMaps, Images
map = applycolourmap(output.frames[1].data, cmap("L4"))
convert.(RGB24, (colorview(RGB, map[:,:,1], map[:,:,2], map[:,:,3])))
