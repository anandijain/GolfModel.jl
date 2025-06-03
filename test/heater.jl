using Setfield
@variables temp(t)
params = @parameters furnace_on_threshold = 0.5 furnace_off_threshold = 0.7 furnace_power = 1.0 leakage = 0.1 furnace_on(t)::Bool = false
eqs = [
    D(temp) ~ furnace_on * furnace_power - temp^2 * leakage
]
furnace_disable = ModelingToolkit.SymbolicContinuousCallback(
    [temp ~ furnace_off_threshold],
    ModelingToolkit.ImperativeAffect(modified=(; furnace_on)) do x, o, c, i
        @set! x.furnace_on = false
    end)
furnace_enable = ModelingToolkit.SymbolicContinuousCallback(
    [temp ~ furnace_on_threshold],
    ModelingToolkit.ImperativeAffect(modified=(; furnace_on)) do x, o, c, i
        @set! x.furnace_on = true
    end)

[temp ~ furnace_off_threshold] => ModelingToolkit.ImperativeAffect(modified=(;
    furnace_on)) do x, o, i, c
    @set! x.furnace_on = false
end

@named sys = ODESystem(
    eqs, t, [temp], params; continuous_events=[furnace_disable, furnace_enable])
ss = structural_simplify(sys)
prob = ODEProblem(ss, [temp => 0.0, furnace_on => true], (0.0, 10.0))
sol = solve(prob)
plot(sol)
hline!([sol.ps[furnace_off_threshold], sol.ps[furnace_on_threshold]],
    l=(:black, 1), primary=false)