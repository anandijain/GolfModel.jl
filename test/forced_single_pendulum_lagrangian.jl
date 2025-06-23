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

continuous_events = [th1 ~ 0] => (stop_affect!, (;))
@named single = lagrangian2system(single_pendulum_lagrangian, qdot, q, p, t, D; Q, continuous_events)

@assert ModelingToolkit.validate(equations(single))
@assert !isempty(ModelingToolkit.continuous_events(single))

# initial conditions and parameters
u0_p_dict = Dict([
    # initial conditions
    th1 => -pi / 2,
    # parameters
    l1 => 1.8 , #meters
    m1 => .4, # using MoI
    tau => 200, #(Lampsa 1975)
    g => 9.80665
])


## Simulation
prob = ODEProblem(single, u0_p_dict, (0.0, 3.0))
sol = solve(prob; saveat=0.005)
plot(sol)

x1dot_expr = l1*cos(th1)*om1
y1dot_expr = l1*sin(th1)*om1
v1_expr = sqrt(x1dot_expr^2 + y1dot_expr^2)
plot(sol.t, sol[v1_expr])
final_club_head_veloctiy = sol[v1_expr][end] # m/s 
anim_single_pend("forced_single_pendulum_lagrangian.mp4",sol,single, :th1)