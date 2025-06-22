using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: t_nounits as tn, D_nounits as Dn
using Plots, Printf

@variables th1(tn) om1(tn) 
@parameters l1 m1 g

x1 = l1 * sin(th1)
y1 = -l1 * cos(th1)

# x1_dot = l1*cos(th1)*om1

function L((om1,), (th1, ), (l1, m1, g), t)

    v1_sq = l1^2 * om1^2

    # kinetic
    K = (1 / 2) * (m1 * v1_sq)

    y1_from_rest = l1 * (1 - cos(th1))

    # potential
    P = g * m1 * y1_from_rest 
    K - P
end
lag_single = L((om1,), (th1, ), (l1, m1, g), t)

function lagrangian2system(
    L, qdot, q, p, t;
    Q=zeros(length(q)),
    defaults=[qdot; q] .=> 0.0,
    kwargs...
)
    Q_vals = Q
    inds = eachindex(q)

    @variables v[inds] x[inds] Q(t)[inds]
    sub = Base.Fix2(substitute, Dict([collect(v .=> qdot); collect(x .=> q)]))
    Lf = L(v, x, p, t)

    F = ModelingToolkit.gradient(Lf, x) + Q
    Lv = ModelingToolkit.gradient(Lf, v)
    rhs = Num.(collect(sub.(F - ModelingToolkit.jacobian(Lv, x) * qdot - ModelingToolkit.derivative.(Lv, (t,)))))
    M = sub.(ModelingToolkit.jacobian(Lv, v))

    eqs = [
        ModelingToolkit.scalarize(D.(qdot) .~ M \ rhs)
        ModelingToolkit.scalarize(D.(q) .~ qdot)
        collect(Q .~ Q_vals)
    ]
    sys = ODESystem(eqs, t, [qdot; q; Q], p; defaults=defaults, kwargs...)
    structural_simplify(sys)
end

# from the original lag 
# # Cart input force
# F = 1000sin(t)

# # Generalized forces
# Q = [F, 0]

# Cart input force
F = 200 # N*m (Lampsa 1975)

# Generalized forces
Q = [F]

q = [th1]
qdot = [om1]
p = [l1, m1, g]

function stop_affect!(integ, u, p, ctx)
    if integ.t <= 0
        return
    else
        terminate!(integ)
    end
end
continuous_events = [th1 ~ 0] => (stop_affect!, [], [], [], nothing)

@named single = lagrangian2system(L, qdot, q, p, tn; Q, continuous_events)
# ModelingToolkit.continuous_events(single)
ic = [
    th1 => -pi / 2,
]

# Parameters
p = Dict([
    l1 => 1.8 , #meters
    m1 => .4, # using MoI
    g => 9.80665
])


## Simulation
prob = ODEProblem(single, ic, (0.0, 3.0), [p...])
sol = solve(prob; saveat=0.005)
plot(sol)

x1dot_expr = l1*cos(th1)*om1
y1dot_expr = l1*sin(th1)*om1
v1_expr = sqrt(x1dot_expr^2 + y1dot_expr^2)
plot(sol.t, sol[v1_expr])
final_club_head_veloctiy = sol[v1_expr][end] # m/s 
anim_single_pend("forced_single_pendulum_lagrangian.mp4",sol,single, :th1)