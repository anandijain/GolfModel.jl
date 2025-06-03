using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Plots

@mtkmodel harm_osc begin
    @parameters begin
        m=1
        k=1
    end
    @variables begin
        x(t) = 1.0 
    end
    @equations begin
        (D^2)(x) ~ -k/m * x
    end
end

@mtkbuild sys = harm_osc()

prob = ODEProblem(sys, [D(sys.x) =>0], (0.0, 10.0), [])
sol = solve(prob)
plot(sol)
plot(sol, idxs=[sys.x, D(sys.x)])

