module GolfModel
using ModelingToolkit, DifferentialEquations, Plots, Printf
using ModelingToolkit: t, D
using DynamicQuantities
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
    System(eqs, t, [qdot; q; Q], p; defaults=defaults, kwargs...)
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
"""
Warning this assumes that v2 is 0
"""
function compute_elastic_collision(m1, v1, m2)
    A = m1 + m2
    B = -2 * m1 * v1
    C = -(m2 - m1) * v1^2

    # we disregard the solution where v1 = v1p, as this means no collision occured. 
    v1p = quadratic_solve(A, B, C)
    v1p_val = v1p[findfirst(x -> !isapprox(x, v1), v1p)]

    # more algebra and we get the golf ball veloctiy FINALLY 
    v2p_val = (m1 / m2) * (v1 - v1p_val)
    v1p_val, v2p_val
end
@mtkmodel simple_ball_model begin
    @description """
    Given an initial velocity and angle, we solve for the trajectory (x(t), y(t)).
    
    I will attempt to end the simulation (t>0) and y(t) = 0.
            """
    @parameters begin
        loft = deg2rad(14), [description = "Initial angle of the ball in radians [rad] from horizontal", unit = u"rad"]
        g = 9.81, [description = "Acceleration due to gravity in meters per second squared [m/s^2]", unit = u"m/s^2"]
        v0 = 76.0, [description = "Initial velocity of the ball in meters per second [m/s]", unit = u"m/s"]
        # k = 0.0052, [description = "drag coefficient of the ball", unit = u"1/m"]
        k = 0, [description = "Drag coefficient of the ball", unit = u"1/m"]
    end
    @variables begin
        x(t) = 0.0, [description = "Horizontal position of the ball in meters [m]", unit = u"m"]
        vx(t) = v0 * cos(loft), [description = "Horizontal velocity of the ball in meters per second [m/s]", unit = u"m/s"]
        y(t) = 0.0, [description = "Vertical position of the ball in meters [m]", unit = u"m"]
        vy(t) = v0 * sin(loft), [description = "Vertical velocity of the ball in meters per second [m/s]", unit = u"m/s"]
        v(t), [description = "Magnitude of the ball's velocity in meters per second [m/s]", unit = u"m/s"]
    end
    @equations begin
        v ~ sqrt(vx^2 + vy^2)
        D(x) ~ vx
        D(vx) ~ -k * v * vx
        D(y) ~ vy
        D(vy) ~ -k * v * vy - g
    end
    @continuous_events begin
        [y ~ 0] => (stop_affect!, (;))
    end

end

function dubble(;kwargs...)
    @variables th1(t) [unit = u"rad", description = "Shoulder rotation from vertical"] om1(t) [unit = u"rad/s", description = "Shoulder rotational velocity"] th2(t) [unit = u"rad", description = "Wrist rotation from vertical"] om2(t) [unit = u"rad/s", description = "Wrist rotational velocity"]
    @parameters l1 [unit = u"m", description = "Arm length"] m1 [unit = u"kg", description = "Arm mass"] l2 [unit = u"m", description = "Club length"] m2 [unit = u"kg", description = "Club mass"] g [unit = u"m/s^2", description = "Gravitational acceleration"]
    @parameters tau_sh [unit = u"N*m", description = "Maximum shoulder torque"] tau_wr [unit = u"N*m", description = "Maximum wrist torque"]
    Symbolics.setmetadata(t, ModelingToolkit.VariableDescription, "Time")

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
        tau_sh => 146.9,
        tau_wr => 15,
        t => 0
        # trel => .15
    ])

    # is there a way in MTK to say rootfind on th1~0 and th2~0- no i dont think so 
    continuous_events = [th1 ~ th2] => (stop_affect!, (;))
    @named dubbl = lagrangian2system(double_pendulum_lagrangian, qdot, q, p, t, D; Q, defaults=defs, continuous_events, kwargs...)
end
function model_syms(sys)
    sts = unknowns(sys)
    ps = parameters(sys)
    ivs = independent_variables(sys)
    [ivs..., sts..., ps...]
end

export lagrangian2system, anim_single_pend, single_pendulum_lagrangian, double_pendulum_lagrangian, stop_affect!
export compute_elastic_collision, model_syms
end # module GolfModel
