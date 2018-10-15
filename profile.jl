# Model profiling

using ProfileView, Profile

Profile.clear()
@profile sim!(output, model, init, layers; time=20)
ProfileView.view()
