using Plots

# your solution arrays
θ1 = sol[th1]
θ2 = sol[th2]
ts = sol.t
l1_val = defs[l1]
l2_val = defs[l2]

# choose every k-th frame to plot, or pick a fixed number of samples
k = 2                      # every 20th timestep
sample_idxs = 1:k:length(ts)

# start a fresh plot
plt = plot(
    xlims=(-(l1_val + l2_val), l1_val + l2_val),
    ylims=(-(l1_val + l2_val), l1_val + l2_val),
    aspect_ratio=:equal,
    legend=false,
    title="Double-pendulum strobe (every $k-th frame)"
)

# overlay each sampled configuration
for (j, i) in enumerate(sample_idxs)
    th1i, th2i = θ1[i], θ2[i]
    x1 = l1_val * sin(th1i)
    y1 = -l1_val * cos(th1i)
    x2 = x1 + l2_val * sin(th2i)
    y2 = y1 - l2_val * cos(th2i)

    # use decreasing alpha so earlier strobes are fainter
    α = 0.1 + 0.9 * (j / length(sample_idxs))

    plot!(
        plt,
        [0, x1, x2],
        [0, y1, y2],
        lw=2,
        alpha=α,
        c=:black,
    )
end

# optionally highlight the very last pose in color
last_i = sample_idxs[end]
x1f = l1_val * sin(θ1[last_i])
y1f = -l1_val * cos(θ1[last_i])
x2f = x1f + l2_val * sin(θ2[last_i])
y2f = y1f - l2_val * cos(θ2[last_i])
plot!(plt, [0, x1f, x2f], [0, y1f, y2f], lw=3, c=:red)
scatter!(plt, [x1f, x2f], [y1f, y2f], ms=8, c=[:blue :red])

# save or display
savefig(plt, "double_pendulum_strobe.png")
