using Plots
plotly()

include("state.jl")
include("simulator.jl")
include("policy.jl")
include("params.jl")

params = Params(0.75, 0.5, 500, 0.05, 1, 5, 5, 1)
initial_state = State(10, 100, 10, 10, 1)

n = 50

sim(initial_state, params, n)
