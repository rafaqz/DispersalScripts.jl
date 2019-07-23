using LinuxPerf
import LinuxPerf: make_bench, enable!, disable!, reset!, reasonable_defaults, counters
const bench = make_bench(reasonable_defaults);

tstop=68
output = ArrayOutput(init, tstop)

@noinline function cacheperf(output, ruleset, tstop)
   enable!(bench)
   sim!(output, ruleset; tstop=tstop)
   disable!(bench)
end

cacheperf(output, ruleset, tstop)
counters(bench)
reset!(bench)
