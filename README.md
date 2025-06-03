# MATH2A (contrary to popular belief)

Golf project plan:

1) single unforced pendulum with point mass 
    - todo: add coupling to ball collision, momentum transfer
    - has_collided flag in system? how do we log it and plot it 
        - "discrete variable" is what they call this 
        - plotting is related to https://github.com/SciML/ModelingToolkit.jl/issues/3010
    - add an effect that on the first collision, we reduce the velocity of the club by the correct amount 
    - figure out how to end the simulation with an event?
    - for the ball model, we cant pick a tspan we need to run it until the affect terminates the simulation
    - todo: fix all of the parameters, but especiial the moment of inertia, as we treat the club as a point mass
2) single forced pendulum 
    - this mimics the shoulder adding a constant force to the swing 
3) double pendulum
