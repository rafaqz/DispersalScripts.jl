# Model profiling

# ProfileView
using ProfileView, Profile, BenchmarkTools
Profile.clear()
@profiler sim!(output, ruleset)
ProfileView.view()


# Juno
tstop=68
output = ArrayOutput(init, tstop)
output = Dispersal.RegionOutput(init, tstop, objective)
data = CellularAutomataBase.simdata(ruleset, init)
@time sim!(output, ruleset; tstop=tstop, data=data)
@time CellularAutomataBase.initdata!(data, ruleset, init)
@time CellularAutomataBase.simloop!(output, ruleset, data, 2:tstop)
@time CellularAutomataBase.swapsource(data)
@time CellularAutomataBase.maprule!(data, ruleset.rules[1])
@time for i in 1:100 sim!(output, ruleset; tstop=tstop, data=data) end
@profiler for i in 1:10 sim!(output, ruleset; tstop=tstop, data=data) end

using Profile, Coverage
Profile.clear_malloc_data()
dirnames = "../CellularAutomataBase/src"
info = analyze_malloc(dirnames)
display(info[end-10:end])
