#=
#Visualise and test status blocks.
framework.jl needs this change around line 118
isshowable(output, t) && showframe(output, ruleset, output[end], t, data)
=#

using Colors
using Gtk
using Cairo
using CellularAutomataGtk
import CellularAutomataBase: showframe
CellularAutomataBase.showframe(o::GtkOutput, ruleset::AbstractRuleset, frame::AbstractArray, t, data) = begin
    # Cairo shows images permuted
    normed = CellularAutomataBase.normaliseframe(ruleset, frame)
    img = similar(normed, RGB24)
    r = CellularAutomataBase.maxradius(ruleset.rules)
    blocksize = 2r
    for i in CartesianIndices(normed)
        blockindex = CellularAutomataBase.indtoblock.((i[1] + r,  i[2] + r), blocksize)
        state = normed[i]
        img[i] = if CellularAutomataBase.sourcestatus(data)[blockindex...]
            if state > 0
                RGB24(state)
            else
                RGB24(0.0, 0.5, 0.5)
            end
        elseif state > 0
            RGB24(1.0, 0.0, 0.0)
        else
            RGB24(0.5, 0.5, 0.0)
        end
    end
    img = permutedims(img)
    println(t)
    Gtk.@guarded Gtk.draw(o.canvas) do widget
        ctx = Gtk.getgc(o.canvas)
        Cairo.image(ctx, Cairo.CairoImageSurface(img), 0, 0,
                    Graphics.width(ctx), Graphics.height(ctx))
    end
end
