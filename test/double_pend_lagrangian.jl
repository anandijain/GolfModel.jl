using GolfModel
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t, D
using DynamicQuantities
using Plots, Printf

@variables th1(t) [unit = u"rad"] om1(t) [unit = u"rad/s"] th2(t) [unit = u"rad"] om2(t) [unit = u"rad/s"]
@parameters l1 [unit = u"m"] m1 [unit = u"kg"] l2 [unit = u"m"] m2 [unit = u"kg"] g [unit = u"m/s^2"]
@parameters tau_sh [unit = u"N*m"] tau_wr [unit = u"N*m"] trel [unit=u"s"]

x1 = l1 * sin(th1)
y1 = -l1 * cos(th1)
x2 = x1 + l2 * sin(th2)
y2 = y1 - l2 * cos(th2)
v1_sq = l1^2 * om1^2
v2_sq = v1_sq + l2^2 * om2^2 + 2 * l1 * l2 * om1 * om2 * cos(th1 - th2)

dpl = double_pendulum_lagrangian((om1, om2), (th1, th2), (l1, m1, l2, m2, g), t)

# tau_wr_t = ifelse(t < trel, 0, tau_wr)
Q = [tau_sh, tau_wr]

q = [th1, th2]
qdot = [om1, om2]
p = [l1, m1, l2, m2, g, tau_sh, tau_wr]

defs = Dict([
    th1 => -pi / 2,
    om1 => 0,
    th2 => -pi,
    om2 => 0,
    l1 => 0.7,
    m1 => 1.6,
    # club params
    l2 => 1.1,
    m2 => 0.34,
    g => 9.80665,
    tau_wr => 30,
    tau_sh => 80,
    # trel => .15
]) 

# is there a way in MTK to say rootfind on th1~0 and th2~0- no i dont think so 
continuous_events = [th1 ~ th2] => (stop_affect!, (;))
@named dubbl = lagrangian2system(double_pendulum_lagrangian, qdot, q, p, t, D; Q, defaults=defs, continuous_events)
@assert ModelingToolkit.validate(equations(dubbl))
@assert !isempty(ModelingToolkit.continuous_events(dubbl))

prob = ODEProblem(dubbl, defs, (0.0, 4.0))
sol = solve(prob; saveat=0.01)
plot(sol)
vel_t = sol[sqrt(v2_sq)]
vel_plot = plot(sol.t, vel_t)

θ1 = sol[th1]
θ2 = sol[th2]
ts = sol.t

l1_val = defs[l1]
l2_val = defs[l2]

anim = @animate for i in eachindex(ts)
    th1i, th2i = θ1[i], θ2[i]

    x1 = l1_val * sin(th1i)
    y1 = -l1_val * cos(th1i)
    x2 = x1 + l2_val * sin(th2i)
    y2 = y1 - l2_val * cos(th2i)

    plot([0, x1, x2], [0, y1, y2],
        lw=2, c=:black, legend=false,
        xlims=(-(l1_val + l2_val), l1_val + l2_val),
        ylims=(-(l1_val + l2_val), l1_val + l2_val),
        aspect_ratio=:equal,
        title=@sprintf("t = %.2f s", ts[i]))
    scatter!([x1, x2], [y1, y2], ms=8, c=[:blue :red])
end

# save out a gif @ 30 fps
mp4(anim, "forced_double_pendulum_ic_90_90_q_80_timed_30.mp4", fps=60)