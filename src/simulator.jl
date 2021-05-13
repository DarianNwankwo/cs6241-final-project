include("state.jl")
include("params.jl")
include("policy.jl")

function step(s::State, p::Params, po::Policy)
    # set up M and F
    M = (1-p.α)s.R + s.MW
    F = p.α*s.R + s.FW

    R = s.R + p.dt * (p.r * (1 - ((F + M) / p.K)) * ((p.α*s.R*(1-p.α)*s.R)/(p.b+M)) - p.δ*s.R)
    FW = s.FW + p.dt * (p.r * p.α * (1 - ((F + M) / p.K)) * s.FW * (M / (p.b + M)) - p.δ * s.FW + po.αWF)
    MW = s.MW + p.dt * (p.r * (1 - p.α) * (1 - ((F + M) / p.K)) * s.FW * (M / (p.b + M)) - p.δ * s.MW + po.αWM)

    # new_state = State(FR, MR, FW, MW, s.t+1)
    new_state = State(R, FW, MW, s.t+1)

    return new_state
end

function sample(is::State, p::Params, n::Int64)
    state = is
    states = State[]
    actions = Policy[]
    push!(states, state)
    for i = 2:n
        αWM = rand() * p.αWMMax
        αWF = rand() * p.αWFMax
        po = Policy(αWM, αWF)
        state = step(state, p, po)
        push!(actions, po)
        push!(states, state)
    end

    return states, actions
end

function sim(is::State, p::Params, po::Policy, n::Int64)
    state = is
    states = State[]
    push!(states, state)
    for i = 2:n
        state = step(state, p, po)
        push!(states, state)
    end

    return states
end

function generate_trajectories(is::State, p::Params, n::Int64, num_trajectories::Int64)
    trajectories = []

    for i = 1:num_trajectories
        push!(trajectories, sample(is, p, n))
    end

    return trajectories
end
