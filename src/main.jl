using Plots
plotly()

include("state.jl")
include("simulator.jl")
include("policy.jl")
include("params.jl")

params = Params(0.75, 0.5, 500, 0.05, 1, 5, 5, 1)
# initial_state = State(10, 100, 10, 10, 1)
initial_state = State(110, 10, 10, 1)

n = 50

states = sim(initial_state, params, n)

tvals = [s.t for s in states]
# frs = [s.R*p.α for s in states]
# mrs = [s.MR*(1-p.α) for s in states]
rs = [s.R for s in states]
fws = [s.FW for s in states]
mws = [s.MW for s in states]

# plot(tvals, frs, label="FR")
# plot!(tvals, mrs, label="MR")
plot(tvals, rs, label="R")
plot!(tvals, fws, label="FW")
plot!(tvals, mws, label="MW")


savefig("testplot.png")
