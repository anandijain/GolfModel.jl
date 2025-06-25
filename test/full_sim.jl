@mtkcompile dubbl = GolfModel.dubble()
@mtkcompile ball = GolfModel.simple_ball_model()

dsyms = model_syms(dubbl)
bsyms  = model_syms(ball)
intersect(dsyms, bsyms)
@unpack th1, om1, th2, om2, l1, m1, l2, m2, g, tau_sh, tau_wr = dubbl
@unpack v0, loft, v, x, y, vx, vy, k = ball

# we want to sim thru the double pendulum and then the ball
model_updates = Dict([tau_sh=>90])
prob1 = ODEProblem(dubbl, model_updates, (0, 5))
sol1 = solve(prob1; saveat=0.01)
v1_sq = l1^2 * om1^2
v2_sq = v1_sq + l2^2 * om2^2 + 2 * l1 * l2 * om1 * om2 * cos(th1 - th2)
clubhead_v = sol1[sqrt(v2_sq)][end]
defs = merge(ModelingToolkit.defaults(dubbl), model_updates)
clubv, ball_v = compute_elastic_collision(defs[m2], clubhead_v, 0.045)

prob2 = ODEProblem(ball, [v0 => ball_v, k => 0.0052], (0, 5))
sol2 = solve(prob2; saveat=0.01)
sol2[x][end]

prob3 = ODEProblem(ball, [v0 => ball_v, k => 0.0052], (0, 5))
sol3 = solve(prob3; saveat=0.01)
sol3[x][end]

traj_plot = plot(sol2[x], sol2[y]; title="Predicted Trajectory for arm_flex_l", yaxis="Height (m)", xaxis="Distance (m)", label="Trajectory", xlabel="Distance (m)", ylabel="Height (m)", legend=:topright)
traj2 = plot!(traj_plot, sol3[x], sol3[y]; label="Trajectory with drag", xlabel="Distance (m)", ylabel="Height (m)")
savefig(traj_plot, "arm_flex_l_trajectory_drag.png")