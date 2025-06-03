@variables x(t)
@parameters c(t)

@mtkbuild sys = ODESystem(
    D(x) ~ c * cos(x), t, [x], [c]; discrete_events=[1.0 => [c ~ c + 1]])

prob = ODEProblem(sys, [x => 0.0], (0.0, 2pi), [c => 1.0])
sol = solve(prob, Tsit5())
plot(sol, idxs=(t, :x, :c))


plot(sol, idxs=:c)
sol.t
plot(sol.t, sol[c])
sol([1.5, 2.5], idxs=[c, c * cos(x)])


using ModelingToolkit, OrdinaryDiffEq, Plots
t = ModelingToolkit.t_nounits
D = ModelingToolkit.D_nounits
@variables x(t)
@parameters c(t)
@mtkbuild sys = ODESystem(
    D(x) ~ c * cos(x), t, [x], [c]; discrete_events=[1.0 => [c ~ c + 1]])
prob = ODEProblem(sys, [x => 0.0], (0.0, 2pi), [c => 1.0])
sol = solve(prob, Tsit5())
plot(sol; idxs=[sys.c, sys.x])