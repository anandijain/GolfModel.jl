using BifurcationKit, Plots
F(x, p) = @. x^2*(x+p[1])
prob = BifurcationProblem(F, [0.], [-1.0], 1;
    record_from_solution=(x, p; k...) -> x[1])
br = continuation(prob, PALC(), ContinuationPar(p_min=-1.0, p_max=1.0))
plot(br)
