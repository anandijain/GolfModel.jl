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
        vx(t) = v0*cos(loft), [description = "Horizontal velocity of the ball in meters per second [m/s]", unit = u"m/s"]
        y(t) = 0.0, [description = "Vertical position of the ball in meters [m]", unit = u"m"]
        vy(t) = v0*sin(loft), [description = "Vertical velocity of the ball in meters per second [m/s]", unit = u"m/s"]
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

@mtkcompile sys_ball = simple_ball_model()

prob_ball = ODEProblem(sys_ball, [sys_ball.v0=>vel_t[end]], (0, 100))
sol = sol_ball = solve(prob_ball; saveat=0.01)
plot(sol_ball[x], sol_ball[y])
# plot(sol_ball[[x,y]])
sol[x][end]