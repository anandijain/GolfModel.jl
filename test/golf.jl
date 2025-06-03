using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using DynamicQuantities
using Plots
using DataFrames, CSV
using Symbolics, Groebner

@show "usings done"
"""
Random notes:

This documentation describes the language for building models. 
https://docs.sciml.ai/ModelingToolkit/stable/basics/MTKLanguage/

good example 
https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/ 
"""


@mtkmodel unforced_single_pendulum begin
    @description  """
            This model treats the mass of the club as a point mass at the club head. 

            In the case of the primitive unforced model, we can algebraically solve for the final club head velocity. 
            F = m * a = m * g
            If θ₀ is the initial angle of the club relative to vertical down, then the change in height of the club head is:
            Δh = L * (1 - cos(θ₀))
            The potential energy at the top is m * g * Δh, and the kinetic energy at the bottom is (1/2) * m * v².
            So the final velocity v at the bottom is given by:
            v = sqrt(2 * g * Δh)) = sqrt(2 * g * L * (1 - cos(θ₀)))

            But this model gets solved numerically, 
            """
    @parameters begin
        L = 1.0, [description = "Length of the golf club in meters [m]", unit = u"m"]
        m = 0.2, [description = "Mass of the golf club in kilograms [kg]", unit = u"kg"]
        g = 9.81, [description = "Acceleration due to gravity in meters per second squared [m/s^2]", unit = u"m/s^2"]

        # gets cancelled out by the assumption that the mass is a point mass at the club head
        # I = m * L^2, [description = "Moment of inertia of the golf club in kg*m^2", unit = u"kg*m^2"]
        has_collided(t)::Bool = false, [description = "Flag to indicate if the club has collided with the ball"]
    end
    @variables begin
        # has_collided(t) = 0, [description = "This is a flag to indicate if the club has collided with the ball"]
        θ(t) = π / 2, [description = "This is the angle of the pendulum from vertical (down) in radians [rad]", unit = u"rad"]
        ω(t) = 0.0, [description = "This is the angular velocity of the pendulum in radians per second [rad/s]", unit = u"rad/s"]
    end
    @equations begin
        # I*(D(θ)^2) ~ -m * g * L * sin(θ) # full thing, but I cancels out to be 
        D(θ) ~ ω # this is the angular velocity
        D(ω) ~ -g / L * sin(θ) # this is the simplified versions
    end
    @continuous_events begin
        # this forces the solver to step when theta is zero, allowing us to get the velocity 
        [θ ~ 0] => [has_collided ~ true]
        [θ ~ 0] => [θ ~ 0]

    end
end

@mtkbuild sys1 = unforced_single_pendulum()
# NOTE this assumed that the initial velocity was zero, which is okay.
prob = ODEProblem(sys1, [], (0.0, 10.0), [])
sol_golf = sol = solve(prob)

# this doesn't really make sense, but it does show that we can index has_collided
# it cant be plotted wrt t easily from what i can tell
sol_golf[sys1.has_collided]

# because of https://github.com/SciML/ModelingToolkit.jl/issues/3010
# you need to plot it with something else for some reason 
plot(sol_golf, idxs=[sys1.θ, sys1.has_collided])

df = DataFrame(sol)
CSV.write("golf.csv", df)
plot(sol)
@show "ode solve done"

# we see an ellipse in theta-omega plane
plot(sol, idxs=(:θ, :ω))
plot(sol, idxs=(t, :θ, :ω))

"""
Now we can take the velocity at the bottom of the swing and use it to simulate the ball. 
"""
zeros = sol(sol.t[sol[sys1.θ==0]])
velocity = abs(zeros[1][2]) # m/s 

# we let u1, u2 be the velocity of the club and ball resp. before collision
# v1,v2 the club and ball velocity after collision
m1 = sys1.m
ModelingToolkit.getdefault(m1)

m1 = 0.2 #club 
m2 = 0.045 #ball
u1 = velocity
@variables v1 v2
# collision eqs

"""
These equations represent conservation of momentum and kinetic energy of the system.

This allows to determine what the velocities of the club and ball are after the collision, given their initial velocities and masses.
"""
cons_moment = m1 * u1 ~ m1 * v1 + m2 * v2
cons_kinetic = m1 * u1^2 ~ m1 * v1^2 + m2 * v2^2
eqs = [
    cons_moment,
    cons_kinetic
]

# algebraic simplification of the equations gives us the following expression for v1 
# Av1^2 + Bv1 + C = 0
A = m1 + m2
B = -2 * m1 * u1
C = -(m2 - m1) * u1^2

# we disregard the solution where u1 = v1, as this means no collision occured. 
quadratic_solve(a, b, c) = ((-b + sqrt(b^2 - 4 * a * c)) / (2 * a), (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
v1_real = quadratic_solve(A, B, C)
v1_val = v1_real[findfirst(x -> !isapprox(x, u1), v1_real)]

# more algebra and we get the golf ball veloctiy FINALLY 
v2_val = (m1 / m2) * (u1 - v1_val)

@mtkmodel simple_ball_model begin
    @description """
    Given an initial velocity and angle, we solve for the trajectory (x(t), y(t)).
    
    I will attempt to end the simulation (t>0) and y(t) = 0.
           """
    @parameters begin
        ball_mass = 0.045 # kg 
        initial_angle = π/4, [description = "Initial angle of the ball in radians [rad] from horizontal", unit = u"rad"]
        g = 9.81, [description = "Acceleration due to gravity in meters per second squared [m/s^2]", unit = u"m/s^2"]
        initial_velocity = v2_val, [description = "Initial velocity of the ball in meters per second [m/s]", unit = u"m/s"]
    end
    @variables begin
        x(t) = 0.0, [description = "Horizontal position of the ball in meters [m]", unit = u"m"]
        y(t) = 0.0, [description = "Vertical position of the ball in meters [m]", unit = u"m"]
    end
    @equations begin
        D(x) ~ v2_expr * cos(initial_angle) # horizontal motion
        D(y) ~ v2_expr * sin(initial_angle) - g * t # vertical motion with gravity
    end
    @continuous_events begin

    end
end

@mtkbuild ball_sys = simple_ball_model()
# NOTE this assumed that the initial velocity was zero, which is okay.
ball_prob = ODEProblem(ball_sys, [], (0.0, Inf), [])
sol_golf = sol = solve(ball_prob)
plot(sol_golf)