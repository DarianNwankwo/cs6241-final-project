println("CS6241 Final Project Baby!")

include("state.jl")
include("simulator.jl")
include("policy.jl")
include("params.jl")

params = Params(0.001, 0.001, 0.001, 0.001, 0.1, 0.1, 0.1, 1, 0.001, 0.001, 0.001, 0.001, 0.001)
base_state = State(1000, 10, 0, 100, 100, 0, 0, 0)
base_policy = Policy(0, 0)

sim(base_state, params, base_policy, 200)
