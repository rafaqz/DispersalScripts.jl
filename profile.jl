# Model profiling

using ProfileView, Profile

output = ArrayOutput(init) 
Profile.clear()
@profile sim!(output, model, init, layers; tstop=100)
ProfileView.view()

@time sim!(output, model, init, layers; tstop=100)
