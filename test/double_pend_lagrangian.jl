using GolfModel
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t, D
using DynamicQuantities
using Plots, Printf

@mtkcompile dubbl = GolfModel.dubble(; continuous_events=[])
@unpack th1, om1, th2, om2, l1, m1, l2, m2, g, tau_sh, tau_wr = dubbl
defs = Dict([
    th1 => -3pi / 4,
    om1 => 0,
    th2 => -pi,
    om2 => 0,
    l1 => 0.7,
    m1 => 1.6,
    # club params
    l2 => 1.1,
    m2 => 0.34,
    g => 9.80665,
    tau_sh => 0,
    tau_wr => 0,
    t => 0
    # trel => .15
])
@assert ModelingToolkit.validate(equations(dubbl))
@assert !isempty(ModelingToolkit.continuous_events(dubbl))

prob = ODEProblem(dubbl, defs, (0.0, 20.0))
sol = solve(prob; saveat=0.01)
plot(sol)
v1_sq = l1^2 * om1^2
v2_sq = v1_sq + l2^2 * om2^2 + 2 * l1 * l2 * om1 * om2 * cos(th1 - th2)
vel_t = sol[sqrt(v2_sq)]
vel_plot = plot(sol.t, vel_t; title="arm_flex_r torque driven double pendulum", xlabel="Time (s)", ylabel="Velocity (m/s)", label="Club head velocity", legend=:topright)
savefig(vel_plot, "arm_flex_r_torque_driven_double_pendulum_velocity.png")

vel_t[end]

prob_l = remake(prob; p=[tau_sh=>90.2])
sol_l = solve(prob_l; saveat=0.01)

# extract the “club-head” speed for each
vel_r = sol_r[sqrt(v2_sq)]   # right arm
vel_l = sol_l[sqrt(v2_sq)]   # left  arm

p = plot(sol_r.t, vel_r;
    label="arm_flex_r",
    title="Velocity Comparison",
    xlabel="Time (s)",
    ylabel="Velocity (m/s)",
    legend=:topright)

plot!(p, sol_l.t, vel_l; label="arm_flex_l")

savefig(p, "velocity_comparison.png")

club_v, ball_v = compute_elastic_collision(double_defs[m2], vel_t[end], 0.045)
prob_ball = ODEProblem(sys_ball, [sys_ball.v0 => ball_v], (0, 100))
sol_ball = solve(prob_ball; saveat=0.01)
plot(sol_ball[x], sol_ball[y])
# plot(sol_ball[[x,y]])
sol_ball[x][end]