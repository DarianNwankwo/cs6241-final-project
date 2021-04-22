include("state.jl")
include("params.jl")
include("policy.jl")

function step(s::State, p::Params)
    # set up M and F
    M = (1-p.α)s.R + s.MW
    F = p.α*s.R + s.FW

    # FR = round(Int64, s.FR + p.dt * (p.r * p.α * (1 - ((F + M) / p.K)) * s.FR * (s.MR / (p.b + M)) - p.δ * s.FR))
    # MR = round(Int64, s.MR + p.dt * (p.r * (1 - p.α) * (1 - ((F + M) / p.K)) * s.FR * (s.MR / (p.b + M)) - p.δ * s.MR))
    R = round(Int64, s.R + p.dt * (p.r * (1 - ((F + M) / p.K)) * ((p.α*s.R*(1-p.α)*s.R)/(p.b+M)) - p.δ*s.R))
    FW = round(Int64, s.FW + p.dt * (p.r * p.α * (1 - ((F + M) / p.K)) * s.FW * (M / (p.b + M)) - p.δ * s.FW + p.AWF))
    MW = round(Int64, s.MW + p.dt * (p.r * (1 - p.α) * (1 - ((F + M) / p.K)) * s.FW * (M / (p.b + M)) - p.δ * s.MW + p.AWM))

    # new_state = State(FR, MR, FW, MW, s.t+1)
    new_state = State(R, FW, MW, s.t+1)

    return new_state
end

function sim(is::State, p::Params, n::Int64)
    state = is
    states = State[]
    push!(states, state)
    for i = 2:n
        state = step(state, p)
        push!(states, state)
    end

    return states
end
