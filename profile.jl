# Model profiling

using ProfileView, Profile, BenchmarkTools

Profile.clear()
@profiler sim!(output, ruleset)
ProfileView.view()

src, dst = CellularAutomataBase.allocatestorage(ruleset, output[end])
sze = size(output[end])
data = CellularAutomataBase.simdata(ruleset, src, dst, sze, 1)
rule = humandisp
f(rule, data, init) = begin
    for x in 1:1000
        for i in 1:400, j in 1:400 CellularAutomataBase.applyrule!(rule, data, one(eltype(init)), (i,j)) end
    end
end
@profile f(rule, data, init)
ProfileView.view()
Profile.print()

tstop=68
output = ArrayOutput(init, tstop)
@time sim!(output, ruleset; tstop=tstop)
@time for i in 1:10 sim!(output, ruleset) end
@profiler for i in 1:10 sim!(output, ruleset) end
