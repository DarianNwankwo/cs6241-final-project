include("state.jl")
include("params.jl")
include("policy.jl")

function susceptibility_dynamics(s::State, p::Params)
    return round(Int64, s.SP - (p.αMR*s.MR + p.αFR*s.FR + p.αMW*s.MW + p.αFW*s.FW)*s.SP)
end

function infected_dynamics(s::State, p::Params)
    return round(Int64, s.IP + (p.αMR*s.MR + p.αFR*s.FR + p.αMW*s.MW + p.αFW*s.FW)*s.SP)
end

function recovered_dynamics(s::State, p::Params)
    return round(Int64, s.RP + p.ρ*s.IP)
end

function regular_female_dynamics(s::State, p::Params)
    return round(Int64, s.FR + p.γMRG*s.FR*s.MR - p.δ*s.FR)
end

function regular_male_dynamics(s::State, p::Params)
    return round(Int64, s.MR+p.γFRM*s.FR*s.MR - p.δ*s.MR)
end

function infected_female_dynamics(s::State, p::Params, po::Policy)
    return round(Int64, s.FW + p.γMWG*s.FW*s.MW
        + p.γMRG*s.FR*s.MR + po.αFW - p.δ*s.FW)
end

function infected_male_dynamics(s::State, p::Params, po::Policy)
    return round(Int64, s.MW + p.γFWM*s.FW*s.MW
        + p.γMRF*s.FW*s.MR + po.αMW - p.δ*s.MW)
end

function step(state::State, params::Params, policy::Policy)
    SP = susceptibility_dynamics(state, params)
    IP = infected_dynamics(state, params)
    RP = recovered_dynamics(state, params)
    FR = regular_female_dynamics(state, params)
    MR = regular_male_dynamics(state, params)
    FW = infected_female_dynamics(state, params, policy)
    MW = infected_male_dynamics(state, params, policy)
    n = state.n + 1
    new_state = State(SP, IP, RP, MR, FR, MW, FW, n)

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
