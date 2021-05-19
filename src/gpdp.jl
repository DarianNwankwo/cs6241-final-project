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

function hyperparam_kernel_vgp(states, V)
    s = vcat([[st.R st.FW st.MW] for st in states]...)

    mean = MeanZero()
    best_kernel = SE(0.0, 0.0)
    best_gp = GP(s', V, mean, best_kernel)
    for ls = 0.0:0.25:2.0
        for σ = 0.0:0.25:2.0
            kernel = SE(ls, σ)

            gp = GP(s', V, mean, kernel)
            if gp.mll > best_gp.mll
                best_gp = gp
                best_kernel = kernel
            end
        end
    end

    return best_kernel, best_gp
end

function hyperparam_kernel_qgp(controls, Q)
    c = vcat([[co.αWF co.αWM] for co in controls]...)

    mean = MeanZero()
    best_kernel = SE(0.0, 0.0)
    best_gp = GP(c', Q, mean, best_kernel)
    for ls = 0.0:0.25:2.0
        for σ = 0.0:0.25:2.0
            kernel = SE(ls, σ)

            gp = GP(c', Q, mean, kernel)
            if gp.mll > best_gp.mll
                best_gp = gp
                best_kernel = kernel
            end
        end
    end

    return best_kernel, best_gp
end

function fit_vgp(states, V, kernel)
    s = vcat([[st.R st.FW st.MW] for st in states]...)

    mean = MeanZero()

    gp = GP(s', V, mean, kernel)
    return gp
end

function fit_qgp(controls, Q, kernel)
    c = vcat([[co.αWF co.αWM] for co in controls]...)

    mean = MeanZero()

    gp = GP(c', Q, mean, kernel)
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

    best_gp1 = GP(s', c[:,1], mean, kernel)
    best_gp2 = GP(s', c[:,2], mean, kernel)

    # for ls = 0.0:0.25:4.0
    #     for σ = 0.0:0.25:4.0
    #         gp1 = GP(s', c[:,1], mean, kernel)
    #
    #         if gp1.mll > best_gp1.mll
    #             best_gp1 = gp1
    #         end
    #
    #         gp2 = GP(s', c[:,2], mean, kernel)
    #         if gp2.mll > best_gp2.mll
    #             best_gp2 = gp2
    #         end
    #     end
    # end
    return best_gp1, best_gp2
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

function sample_cumulative_controls(trajectories, n)
    controls = vcat([traj[2][n:end] for traj in trajectories]...)

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

function sample_cumulative_states(trajectories, n)
    states = vcat([traj[1][n:end] for traj in trajectories]...)

    return states
end

function find_min_Q_samples(controls, params, Qgp)
    umin = controls[argmin([predict_y(Qgp, [co.αWF co.αWM]')[1][1] for co in controls])]

    return umin
end

function find_min_Q_grid(controls, params, Qgp)
    grid_controls = [Policy(αWF, αWM) for αWM = 0.0:params.αWMMax/20.0:params.αWMMax for αWF = 0.0:params.αWFMax/20.0:params.αWFMax]

    umin = grid_controls[argmin([predict_y(Qgp, [co.αWF co.αWM]')[1][1] for co in grid_controls])]
    return umin
end

function learn(trajectories, params, n)
    sample_states = params.sample_state_func(trajectories, n)
    V = [terminal_cost(state, params) for state in sample_states]
    # Vkernel, Vgp = hyperparam_kernel_vgp(sample_states, V)
    Vkernel = SE(2.0, 2.0)
    Qkernel = SE(2.0, 2.0)
    Vgp = fit_vgp(sample_states, V, Vkernel)

    # Find best Q kernel from a sample
    # sample_controls = params.sample_control_func(trajectories, n-1)
    # st = sample(sample_states)
    # Qx = [cost(st, co, params) + expectation(Vgp, st, co, params) for co in sample_controls]

    # Qkernel, Qgp = hyperparam_kernel_qgp(sample_controls, Qx)

    final_pis = []
    all_states = []
    push!(all_states, sample_states)
    for i = n-1:-1:1
        pis = []
        V = []

        sample_states = params.sample_state_func(trajectories, i)
        push!(all_states, sample_states)

        # for st in states
        for st in sample_states
            sample_controls = params.sample_control_func(trajectories, i)
            Qx = [cost(st, co, params) + expectation(Vgp, st, co, params) for co in sample_controls]

            Qgp = fit_qgp(sample_controls, Qx, Qkernel)

            umin = find_min_Q_samples(sample_controls, params, Qgp)
            push!(pis, umin)

            cumin = [umin.αWF umin.αWM]'

            push!(V, predict_y(Qgp, cumin)[1][1])
        end

        Vgp = fit_vgp(sample_states, V, Vkernel)
        push!(final_pis, pis)
        # push!(final_pis, learn_pi(sample_states, pis))
        # final_pis = pis
    end

    # return reverse(final_pis)
    return Vgp, reverse(all_states), reverse(final_pis)
end
