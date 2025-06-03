using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Plots

@mtkmodel lotka begin
    @parameters begin
        a = rand()
        b = rand()
        c = rand()
        d = rand()
    end
    @variables begin
        x(t) = 1.0
        y(t) = 1.0
    end
    @equations begin
        D(x) ~ a * x - b * x * y
        D(y) ~ c * x * y - d * y
    end
end

@mtkbuild sys = lotka()

prob = ODEProblem(sys, [], (0.0, 10.0), [])
sol = solve(prob)
plot(sol, idxs=[:x, :y])
plot(sol, idxs=(:x, :y))