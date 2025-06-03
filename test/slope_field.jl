f(t, y) = y - t
f(xs) = xs[2] - xs[1]
vals = 0:10
pts = collect(Iterators.product(vals, vals))
dys = f.(pts)
collect(zip(pts, dys))

function f(du, u, p, t)
    du[1] = 1-cos(u[1]) + (1+cos(u[1]))*p[1]
end

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Plots
prob = ODEProblem(f, [0.0], (0.0, 10.0), [-.1])
plot(solve(prob))
prob2 = remake(prob; p=[0])
plot!(solve(prob2), label="p=0")

foo(x;p) = 1 - cos(x) + (1 + cos(x)) * p
foo.(-5:5;p=0)

f(u, p, t) = 1 - cos(u) + (1 + cos(u)) * p[1]

ts = range(-20, 20, length=20)      # t from 0 to 10
us = range(-2π, 2π, length=20)      # u from -π to π
p = -.1

# 3) Compute slope = du/dt at each grid point
S = [f(u, p, t) for t in ts, u in us]  # S[i,j] = slope at (t=ts[i], u=us[j])
norms = sqrt.(1 .+ S .^ 2)
dx = 0.4 * (1 ./ norms)          # scale each arrow horizontally by 0.4
dy = 0.4 * (S ./ norms)          # scale each arrow vertically by 0.4
quiver(
    repeat(ts, inner=length(us)),      # x‑coords of arrows
    repeat(us, outer=length(ts)),      # y‑coords of arrows
    quiver=(vec(dx), vec(dy)),       # arrow components
    aspect_ratio=1,
    xlabel="t",
    ylabel="u",
    title="Slope field for du/dt = 1 - cos(u) + (1+cos(u))·$(p)"
)

function euler(f, u0, tspan, dt)
    t = tspan[1]
    ts = [t]
    us = [u0]
    u = u0
    while t < tspan[2]
        u += dt * f(u, p, t)
        t += dt
        push!(ts, t)
        push!(us, u)
    end
    (ts, us)
end
f(y, _, _) = 2y
dt = 0.001
ts, xs = euler(f, 1.0, (0.0, 10.0), dt)
idx = Int(.5/dt)
xs[idx]

g(y, _, _) = exp(2/y)
dt = .5
ts, xs = euler(g, 2.0, (0.0, 100), dt)
collect(zip(ts, xs))