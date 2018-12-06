# Model profiling

using ProfileView, Profile

output = ArrayOutput(init, 100) 
Profile.clear()
@profile sim!(output, model, init, layers; tstop=100)
ProfileView.view()

@time sim!(output, model, init, layers; tstop=100)
