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

function learn_pi(states, pis)
    s = vcat([[st.R st.FW st.MW] for st in states]...)
    c = vcat([[co.αWF co.αWM] for co in pis]...)

    mean = MeanZero()
    kernel = SE(2.0, 2.0)

    gp1 = GP(s', c[:,1], mean, kernel)
    gp2 = GP(s', c[:,2], mean, kernel)
    return gp1, gp2
end

function sample_all_controls(trajectories, n)
    controls = vcat([traj[2] for traj in trajectories]...)

    return controls
end

function sample_random_controls(trajectories, n)
    controls = sample_all_controls(trajectories, n)

    return [sample(controls) for j in 1:100]
end

function sample_time_slice_controls(trajectories, n)
    controls = [traj[2][n] for traj in trajectories]

    return controls
end

function sample_all_states(trajectories, n)
    states = vcat([traj[1] for traj in trajectories]...)

    return states
end

function sample_random_states(trajectories, n)
    states = sample_all_states(trajectories, n)

    return [sample(states) for i in 1:1000]
end

function sample_time_slice_states(trajectories, n)
    states = [traj[1][n] for traj in trajectories]

    return states
end

function learn(trajectories, params, n)
    sample_states = params.sample_state_func(trajectories, n)
    V = [terminal_cost(state, params) for state in sample_states]
    Vgp = fit_vgp(sample_states, V)
    final_pis = []
    for i = n-1:-1:1
        pis = []
        V = []

        sample_states = params.sample_state_func(trajectories, i)

        # for st in states
        for st in sample_states
            sample_controls = params.sample_control_func(trajectories, i)
            Qx = [cost(st, co, params) + 0.99*expectation(Vgp, st, co, params) for co in sample_controls]

            Qgp = fit_qgp(sample_controls, Qx)

            #TODO: Numerical methods for finding min instead of looking at support policies
            umin = sample_controls[argmin([predict_y(Qgp, [co.αWF co.αWM]')[1][1] for co in sample_controls])]
            push!(pis, umin)

            cumin = [umin.αWF umin.αWM]'

            push!(V, predict_y(Qgp, cumin)[1][1])
        end

        Vgp = fit_vgp(sample_states, V)
        final_pis = pis
    end

    return Vgp, final_pis, sample_states
end
