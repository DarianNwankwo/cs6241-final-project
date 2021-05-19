using Plots
using GaussianProcesses
using StatsBase
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

function run(n, num_trajectories, name, sample_state_func, sample_control_func)
    params = Params(0.65, 0.25, 500, 0.0025, 1, 10, 10, 0.5, 1, sample_state_func, sample_control_func)
    initial_state = State(30, 0, 0, 1)

    trajectories = generate_trajectories(initial_state, params, n, num_trajectories)

    Vgp, states, pis = learn(trajectories, params, n)

    mkpath(name)

    open(string(name, "/pi_output.txt"), "w") do io
        for p in pis
            write(io, "[\n")
            for co in p
                write(io, string(co.αWF, " ", co.αWM, "\n"))
            end
            write(io, "]\n")
        end
    end

    open(string(name, "/state_output.txt"), "w") do io
        for s in states
            write(io, "[\n")
            for st in s
                write(io, string(st.R, " ", st.FW, " ", st.MW, " ", st.t, "\n"))
            end
            write(io, "]\n")
        end
    end

    pisgp = []
    for i = 1:n-1
        push!(pisgp, learn_pi(states[i], pis[i]))
    end

    states, controls = simgp(initial_state, params, pisgp, n)

    open(string(name, "/opt_trajectories.txt"), "w") do io
        for st in states
            write(io, string(st.FW, " "))
        end
        write(io, "\n")

        for st in states
            write(io, string(st.MW, " "))
        end
        write(io, "\n")

        for st in states
            write(io, string(st.R, " "))
        end
        write(io, "\n")

        for st in states
            write(io, string(st.t, " "))
        end
        write(io, "\n")

        for co in controls
            write(io, string(co.αWF, " "))
        end
        write(io, "\n")

        for co in controls
            write(io, string(co.αWM, " "))
        end
    end

    tvals = [s.t for s in states]
    rs = [s.R for s in states]
    fws = [s.FW for s in states]
    mws = [s.MW for s in states]

    plot(tvals, rs, label="R")
    plot!(tvals, fws, label="FW")
    plot!(tvals, mws, label="MW")

    savefig(string(name, "/plot.png"))

end


run(50, 100, "random_sample", sample_random_states, sample_random_controls)
run(50, 10, "all_sample", sample_all_states, sample_all_controls)
run(50, 100, "time_slice", sample_time_slice_states, sample_time_slice_controls)
run(50, 20, "cumulative_sample", sample_cumulative_states, sample_cumulative_controls)
run(50, 30, "cum_state_rand_controls", sample_cumulative_states, sample_random_controls)
