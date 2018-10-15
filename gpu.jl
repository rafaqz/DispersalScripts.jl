
randomstate = GPUArrays.cached_state(init)
sim!(output, model, init, layers, randomstate; time=10)
