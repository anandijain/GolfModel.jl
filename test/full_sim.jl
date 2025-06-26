using GolfModel
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t, D
using DynamicQuantities
using Plots, Printf
using DataFrames, CSV, SummaryTables


@mtkcompile dubbl = GolfModel.dubble()
@mtkcompile ball = GolfModel.simple_ball_model()

dsyms = model_syms(dubbl)
bsyms = model_syms(ball)
intersect(dsyms, bsyms)
@unpack th1, om1, th2, om2, l1, m1, l2, m2, g, tau_sh, tau_wr = dubbl
@unpack v0, loft, v, x, y, vx, vy, k = ball
v1_sq = l1^2 * om1^2
v2_sq = v1_sq + l2^2 * om2^2 + 2 * l1 * l2 * om1 * om2 * cos(th1 - th2)

nms = ["Amateur", "Intermediate", "Professional"]
sets = [
    [l1 => 0.63, l2 => 1.14, m2 => 0.32, tau_sh => 30, tau_wr => 10],
    [l1 => 0.64, l2 => 1.14, m2 => 0.32, tau_sh => 40, tau_wr => 15],
    [l1 => 0.65, l2 => 1.14, m2 => 0.32, tau_sh => 100, tau_wr => 35]
]
drags = [0, 0.0052]
get_carry(bsol) = bsol[x][end]
get_chvel(ssol) = ssol[sqrt(v2_sq)][end]
sim_names = []
ssols = []
bsols = []
bvels = []

for (i, set) in enumerate(sets)
    for (j, kdrag) in enumerate(drags)
        sim_name = nms[i] * (j == 1 ? " (No Drag)" : " (With Drag)")
        push!(sim_names, sim_name)
        model_updates = Dict(set)

        # we want to sim thru the double pendulum and then the ball
        prob1 = ODEProblem(dubbl, model_updates, (0, 1000)) # sim to 1000 since we have our termination condition
        ssol = solve(prob1; saveat=0.01)
        push!(ssols, ssol)

        clubhead_v = ssol[sqrt(v2_sq)][end]

        defs = merge(ModelingToolkit.defaults(dubbl), model_updates)
        clubv, ball_v = compute_elastic_collision(defs[m2], clubhead_v, 0.045)
        push!(bvels, ball_v)

        prob2 = ODEProblem(ball, [v0 => ball_v, k => kdrag, loft => deg2rad(10)], (0, 1000))
        bsol = solve(prob2; saveat=0.01)
        push!(bsols, bsol)
    end
end
sim_df = DataFrame("Simulation Name" => sim_names, "Clubhead Velocity (m/s)" => get_chvel.(ssols), "Ball Velocity (m/s)" => bvels, "Carry Distance (m)" => get_carry.(bsols))
sim_tbl = simple_table(sim_df)
show(stdout, MIME"text/typst"(), sim_tbl)
traj_plot = plot()
for (i, (ssol, bsol)) in enumerate(zip(ssols, bsols))
    traj_plot = plot!(traj_plot, sol2[x], sol2[y]; title="Predicted Trajectory for arm_flex_l", yaxis="Height (m)", xaxis="Distance (m)", label="Trajectory", xlabel="Distance (m)", ylabel="Height (m)", legend=:topright)
    @printf("Set %d (%s): Clubhead velocity = %.2f m/s\n", i, nms[i], ssol[sqrt(v2_sq)][end])
end
traj2 = plot!(traj_plot, sol3[x], sol3[y]; label="Trajectory with drag", xlabel="Distance (m)", ylabel="Height (m)")
savefig(traj_plot, "arm_flex_l_trajectory_drag.png")

# separate traj validation 
ball_speeds = [59.46, 67.06, 76.25]
bssols = []
dists = []
for ball_v in ball_speeds
    prob2 = ODEProblem(ball, [v0 => ball_v, k => 0], (0, 10))
    bsol = solve(prob2; saveat=0.01)
    push!(bssols, bsol)
    push!(dists, deepcopy(bsol[x][end]))
end
dists