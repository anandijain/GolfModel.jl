module GolfModel
using ModelingToolkit, DifferentialEquations, Plots, Printf

function lagrangian2system(
    L, qdot, q, p, t, D;
    Q=zeros(length(q)),
    defaults=[qdot; q] .=> 0.0,
    kwargs...
)
    Q_vals = Q
    inds = eachindex(q)
    units_x = ModelingToolkit.get_unit.(q)
    units_v = ModelingToolkit.get_unit.(qdot)
    units_Q = ModelingToolkit.get_unit.(Q_vals)

    @variables x[inds] [unit = units_x]
    @variables v[inds] [unit = units_v]
    @variables Q(t)[inds] [unit = units_Q]

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
    sys = System(eqs, t, [qdot; q; Q], p; defaults=defaults, kwargs...)
    mtkcompile(sys)
end

function anim_single_pend(name, sol, sys, theta_sym)
    ts = range(sol.prob.tspan...; step=0.1)
    θ = getproperty(sys, theta_sym)
    s = sol(ts)
    θs = s[θ]      # Array of θ(t) values
    L_val = 1

    # build the animation
    anim = @animate for (i, t) in enumerate(ts)
        θi = θs[i]
        # bob position
        x = L_val * sin(θi)
        y = -L_val * cos(θi)

        # draw the rod and the bob
        plot([0, x], [0, y],
            lw=2, legend=false,
            xlims=(-L_val, L_val),
            ylims=(-L_val, L_val),
            aspect_ratio=:equal,
            title=@sprintf("t = %.2f s", t))
        scatter!([x], [y], ms=8, c=:blue)
    end

    # save out at 30 fps
    mp4(anim, name, fps=60)
end

function single_pendulum_lagrangian((om1,), (th1,), (l1, m1, g), t)

    v1_sq = l1^2 * om1^2

    # kinetic
    K = (1 / 2) * (m1 * v1_sq)

    y1_from_rest = l1 * (1 - cos(th1))

    # potential
    P = g * m1 * y1_from_rest
    K - P
end

function double_pendulum_lagrangian((om1, om2), (th1, th2), (l1, m1, l2, m2, g), t)
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

"""
A helper for continuous_events that terminates the integration on a condition, and enforces that t > 0.
"""
function stop_affect!(mod, obs, ctx, integ)
    if integ.t <= 0
        return (;)
    else
        terminate!(integ)
        (;)
    end
end

quadratic_solve(a, b, c) = ((-b + sqrt(b^2 - 4 * a * c)) / (2 * a), (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))

export lagrangian2system, anim_single_pend, single_pendulum_lagrangian, double_pendulum_lagrangian, stop_affect!

end # module GolfModel
