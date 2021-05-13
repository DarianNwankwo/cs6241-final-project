include("state.jl")
include("simulator.jl")
include("policy.jl")
include("params.jl")

function cost(state, policy, params)
    return (policy.αWF + policy.αWM + state.R * params.α * params.bp) * params.dt
end

function terminal_cost(state, params)
    return 0.0
end

function expectation(gp, state, policy, params)
    ns = step(state, params, policy)

    return predict_y(gp, [ns.R ns.FW ns.MW]')[1][1]
end

function fit_vgp(states, V)
    s = vcat([[st.R st.FW st.MW] for st in states]...)

    mean = MeanZero()
    kernel = SE(2.0, 0.0) # TODO: DON'T HARDCODE

    gp = GP(s', V, mean, kernel)
    return gp
end

function fit_qgp(policies, Q)
    p = vcat([[po.αWF po.αWM] for po in policies]...)

    mean = MeanZero()
    kernel = SE(2.0, 2.0) # TODO: DON'T HARDCODE

    gp = GP(p', Q, mean, kernel)
    return gp
end

function extract_support_states(trajectories)
    states = vcat([traj[1] for traj in trajectories]...)

    return states
end

function extract_support_controls(trajectories)
    controls = vcat([traj[2] for traj in trajectories]...)

    return controls
end

function learn(trajectories, params, n)
    # s = [traj[1][n] for traj in trajectories]
    states = extract_support_states(trajectories)
    controls = extract_support_controls(trajectories)
    V = [terminal_cost(state, params) for state in states]
    Vgp = fit_vgp(states, V)
    final_pis = []
    for i = n-1:-1:1
        pis = []
        V = []
        # s = [traj[1][i] for traj in trajectories]
        # p = [traj[2][i] for traj in trajectories]

        for st in states
            Qx = [cost(st, co, params) + 0.99*expectation(Vgp, st, co, params) for co in controls]

            Qgp = fit_qgp(controls, Qx)

            #TODO: Numerical methods for finding min instead of looking at support policies
            umin = controls[argmin([predict_y(Qgp, [co.αWF co.αWM]')[1][1] for co in controls])]
            push!(pis, umin)

            cumin = [umin.αWF umin.αWM]'

            push!(V, predict_y(Qgp, cumin)[1][1])
        end

        Vgp = fit_vgp(states, V)
        final_pis = pis
    end

    return Vgp, final_pis, states
end
