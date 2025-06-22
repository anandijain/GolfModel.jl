using Printf
# pull out the angle and time vectors
function anim_single_pend(name, sol, sys, theta_sym)
    ts = range(sol.prob.tspan...;step=0.1)
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

anim_single_pend("single_forced_lagrangian.mp4", sol, single, :th1)