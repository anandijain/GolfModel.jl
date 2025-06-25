θ1 = sol[th1]
θ2 = sol[th2]
ts = sol.t

l1_val = defs[l1]
l2_val = defs[l2]

anim = @animate for i in eachindex(ts)
    th1i, th2i = θ1[i], θ2[i]

    x1 = l1_val * sin(th1i)
    y1 = -l1_val * cos(th1i)
    x2 = x1 + l2_val * sin(th2i)
    y2 = y1 - l2_val * cos(th2i)

    plot([0, x1, x2], [0, y1, y2],
        lw=2, c=:black, legend=false,
        xlims=(-(l1_val + l2_val), l1_val + l2_val),
        ylims=(-(l1_val + l2_val), l1_val + l2_val),
        aspect_ratio=:equal,
        title=@sprintf("t = %.2f s", ts[i]))
    scatter!([x1, x2], [y1, y2], ms=8, c=[:blue :red])
end

# save out a gif @ 30 fps
gif(anim, "forced_double_pendulum_ic_HUY.gif", fps=30)