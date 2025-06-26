using GolfModel
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t, D
using DynamicQuantities
using Plots, Printf

@variables th1(t) [unit = u"rad"] om1(t) [unit=u"rad/s"]
@parameters l1 [unit = u"m"] m1 [unit = u"kg"] g [unit = u"m/s^2"]

x1 = l1 * sin(th1)
y1 = -l1 * cos(th1)

lag_single = single_pendulum_lagrangian((om1,), (th1,), (l1, m1, g), t)

# from the original lag 
# F = 1000sin(t)
# Q = [F, 0]

@parameters tau [unit = u"N*m"]
Q = [tau]

q = [th1]
qdot = [om1]
p = [l1, m1, g, tau]
# initial conditions and parameters
defs = Dict([
    # initial conditions
    th1 => -pi / 2,
    om1 => 0,
    # parameters
    l1 => 1.8, #meters
    m1 => 0.4, # using MoI
    tau => 300, #(Lampsa 1975)
    g => 9.80665
])

continuous_events = [th1 ~ 0] => (stop_affect!, (;))
@mtkcompile single = lagrangian2system(single_pendulum_lagrangian, qdot, q, p, t, D; Q, defaults=defs, continuous_events)
@assert ModelingToolkit.validate(equations(single))
@assert !isempty(ModelingToolkit.continuous_events(single))

## Simulation
prob = ODEProblem(single, defs, (0.0, 3.0))
sol = solve(prob; saveat=0.01)
plot(sol, title="Single pendulum driven with 300Nm of Constant Torque")
p = plot!(sol.t, sol[v1_expr], label="Clubhead Velocity (m/s)")

savefig(p, "300nm_single_pendulum.png")

# x1dot_expr = l1*cos(th1)*om1
# y1dot_expr = l1*sin(th1)*om1
# v1_expr = sqrt(x1dot_expr^2 + y1dot_expr^2)
# plot(sol.t, sol[v1_expr])
# final_club_head_veloctiy = sol[v1_expr][end] # m/s 
# anim_single_pend("forced_single_pendulum_lagrangian.mp4",sol,single, :th1)


using Plots
using LaTeXStrings
plot(sol.t, sol[th1], lw=2, label=L"\theta_1\ (rad)")
plot!(sol.t, sol[om1], lw=2, label=L"\omega_1\ (rad/s)")
plot!(sol.t, sol[v1_expr], lw=2, label=L"v_{\rm club}\ (m/s)")

xlabel!("Time (s)")
ylabel!("Quantity")
title!("300 Nm of constant torque - Single pendulum")

# bump up the top margin so the title fits
plot!(margin=Plots.Margins(5mm, 5mm, 15mm, 5mm))
#                       left, right, top, bottom

savefig("pendulum.png")