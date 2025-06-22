using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: t_nounits as tn, D_nounits as Dn
using Plots, Printf

@variables th1(tn) th2(tn) om1(tn) om2(tn)
@parameters l1 m1 l2 m2 g
# D = Differential(t)

x1 = l1 * sin(th1)
y1 = -l1 * cos(th1)
x2 = x1 + l2 * sin(th2)
y2 = y1 - l2 * cos(th2)
v1_sq = l1^2 * om1^2
v2_sq = v1_sq + l2^2 * om2^2 + 2 * l1 * l2 * om1 * om2 * cos(th1 - th2)

# x1_dot = l1*cos(th1)*om1

function L((om1, om2), (th1, th2), (l1, m1, l2, m2, g), t)

    v1_sq = l1^2 * om1^2
    v2_sq = v1_sq + l2^2 * om2^2 + 2 * l1 * l2 * om1 * om2 * cos(th1 - th2)

    # kinetic
    K = (1 / 2) * (m1 * v1_sq + m2 * v2_sq)

    y1_from_rest = l1 * (1 - cos(th1))
    y2_from_rest = y1_from_rest + l2 * (1 - cos(th2))

    # potential
    P = g * (m1 * y1_from_rest + m2 * y2_from_rest)
    K - P
end
L((om1, om2), (th1, th2), (l1, m1, l2, m2, g), t)

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

    Dn = Differential(t)

    eqs = [
        ModelingToolkit.scalarize(D.(qdot) .~ M \ rhs)
        ModelingToolkit.scalarize(D.(q) .~ qdot)
        collect(Q .~ Q_vals)
    ]
    sys = ODESystem(eqs, t, [qdot; q; Q], p; defaults=defaults, kwargs...)
    structural_simplify(sys)
end

T_rel = .15
tau_wrist = 30
tau_wr = ifelse(t < T_rel, 0, tau_wrist)
# F = 0
# Generalized forces
Q = [80, tau_wr]

q = [th1, th2]
qdot = [om1, om2]
p = [l1, m1, l2, m2, g]
function stop_affect!(integ, u, p, ctx)
    if integ.t <= 0
        return
    else
        terminate!(integ)
    end
end
# is there a way in MTK to say rootfind on th1~0 and th2~0
continuous_events = [th2 ~ th1] => (stop_affect!, [], [], [], nothing)
# continuous_events=[]
# Make equations of motion
@named dubbl = lagrangian2system(L, qdot, q, p, tn; Q, continuous_events)

# Initial Conditions
ic = [
    th1 => -pi/2,
    th2 => -pi/2,
]

# Parameters
p = Dict([
    l1 => 0.7,
    m1 => 1.6,
    l2 => 1.1,
    m2 => 0.34,
    g => 9.80665
])


## Simulation
prob = ODEProblem(dubbl, ic, (0.0, 4.0), [p...])
sol = solve(prob; saveat=0.01)
plot(sol)
vel_plot = plot(sol.t, sol[sqrt(v2_sq)])

# extract time series and angles from your solution
θ1 = sol[th1]      # or sol(th1) depending on your DTK version
θ2 = sol[th2]
ts = sol.t

# pull out numeric parameters (in this example we used l1=1, l2=1)
l1_val = p[l1]
l2_val = p[l2]

# build the animation
anim = @animate for i in eachindex(ts)
    # current angles
    th1i, th2i = θ1[i], θ2[i]

    # pendulum coords
    x1 = l1_val * sin(th1i)
    y1 = -l1_val * cos(th1i)
    x2 = x1 + l2_val * sin(th2i)
    y2 = y1 - l2_val * cos(th2i)

    # draw rods + bobs
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