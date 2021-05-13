using Plots
using GaussianProcesses
plotly()

include("state.jl")
include("simulator.jl")
include("policy.jl")
include("params.jl")
include("gpdp.jl")

# include("src/state.jl")
# include("src/simulator.jl")
# include("src/policy.jl")
# include("src/params.jl")
# include("src/gpdp.jl")

params = Params(0.65, 0.25, 500, 0.0025, 1, 10, 10, 0.5, 1)
# initial_state = State(10, 100, 10, 10, 1)
initial_state = State(30, 0, 0, 1)
# policy = Policy(5, 5)

n = 50
num_trajectories = 10

# states = sim(initial_state, params, policy, n)
trajectories = generate_trajectories(initial_state, params, n, num_trajectories)

learn(trajectories, params, n)

# states = trajectories[1][1]
#
# tvals = [s.t for s in states]
# # frs = [s.R*p.α for s in states]
# # mrs = [s.MR*(1-p.α) for s in states]
# rs = [s.R for s in states]
# fws = [s.FW for s in states]
# mws = [s.MW for s in states]
#
# # plot(tvals, frs, label="FR")
# # plot!(tvals, mrs, label="MR")
# plot(tvals, rs, label="R")
# plot!(tvals, fws, label="FW")
# plot!(tvals, mws, label="MW")
#
# savefig("testplot.png")
