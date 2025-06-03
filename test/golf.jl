using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using DynamicQuantities
using Plots

"""
Random notes:

This documentation describes the language for building models. 
https://docs.sciml.ai/ModelingToolkit/stable/basics/MTKLanguage/

good example 
https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/ 
"""

"""
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
@mtkmodel unforced_single_pendulum begin
    @parameters begin
        L = 1.0, [description = "Length of the golf club in meters [m]", unit = u"m"]
        m = 0.2, [description = "Mass of the golf club in kilograms [kg]", unit = u"kg"]
        g = 9.81, [description = "Acceleration due to gravity in meters per second squared [m/s^2]", unit = u"m/s^2"]
        I = m * L^2, [description = "Moment of inertia of the golf club in kg*m^2", unit = u"kg*m^2"]
    end
    @variables begin
        θ(t) = π / 2, [description = "This is the angle of the pendulum from vertical (down) in radians [rad]", unit = u"rad"]
        ω(t) = 0.0, [description = "This is the angular velocity of the pendulum in radians per second [rad/s]", unit = u"rad/s"]
    end
    @equations begin
        # I*(D(θ)^2) ~ -m * g * L * sin(θ) # full thing, but I cancels out to be 
        D(θ) ~ ω # this is the angular velocity
        D(ω) ~ -g / L * sin(θ) # this is the simplified versions
    end
    @continuous_events begin
        [θ ~ 0]
    end
end

@mtkbuild sys1 = unforced_single_pendulum()
# NOTE this assumed that the initial velocity was zero, which is okay.
prob = ODEProblem(sys1, [], (0.0, 10.0), [])
sol = solve(prob)

plot(sol)

# we see an ellipse in theta-omega plane
plot(sol, idxs=(:θ, :ω))
plot(sol, idxs=(t, :θ, :ω))

"""
Now we can take the velocity at the bottom of the swing and use it to simulate the ball. 
"""
sol()
