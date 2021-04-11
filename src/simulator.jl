include("state.jl")
include("params.jl")
include("policy.jl")

function susceptibility_dynamics(state::State, params::Params)
    return round(Int64, state.SP - params.Β*state.SP*state.IMos)
end

function infected_mosquito_dynamics(state::State, params::Params)
    return round(Int64, state.IMos + params.αm*state.MR + params.αf*state.FR)
end

function regular_female_dynamics(state::State, params::Params)
    return round(Int64, state.FR + params.γMRG*state.FR*state.MR)
end

function regular_male_dynamics(state::State, params::Params)
    return round(Int64, state.MR+params.γFRM*state.FR*state.MR)
end

function infected_female_dynamics(state::State, params::Params, policy::Policy)
    return round(Int64, state.FW + params.γMWG*state.FW*state.MW
        + params.γMRG*state.FR*state.MR + policy.αFW)
end

function infected_male_dynamics(state::State, params::Params, policy::Policy)
    return round(Int64, state.MW + params.γFWM*state.FW*state.MW
        + params.γMRF*state.FW*state.MR + policy.αMW)
end

function step(state::State, params::Params, policy::Policy)
    SP = susceptibility_dynamics(state, params)
    IMos = infected_mosquito_dynamics(state, params)
    FR = regular_female_dynamics(state, params)
    MR = regular_male_dynamics(state, params)
    FW = infected_female_dynamics(state, params, policy)
    MW = infected_male_dynamics(state, params, policy)
    n = state.n + 1
    new_state = State(SP, IMos, FR, MR, FW, MW, n)

    return new_state
end

function sim(base_state::State, params::Params, policy::Policy, steps::Int64)
    state = base_state
    println(state)
    for t = 1:steps
        state = step(state, params, policy)
        println(state)
    end
end
