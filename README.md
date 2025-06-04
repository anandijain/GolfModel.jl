# MATH2A (contrary to popular belief)

Golf project plan:

1) single unforced pendulum with point mass 
2) single forced pendulum 
3) double unforced pendulum
4) double forced pendulum 


## Todo: 
- forcing functions https://chatgpt.com/c/683fb440-f8a0-8009-ad99-a5231d7343c5
- coefficient of restitution for the collision
- compose both models heirarchally
- (do after composing models) add an effect that on the first collision, we reduce the velocity of the club by the correct amount 
- fix model parameter values and IC to be realistic with actual golf, particularly club moment of inertia
- units validation 

# Tier 2 todos:
- use the MTKStdLib 
- figure out latexify the eqs 

## Done:
- add coupling to ball collision, momentum transfer
- has_collided flag in system? how do we log it and plot it 
    * "discrete variable" is what they call this 
    * plotting is related to https://github.com/SciML/ModelingToolkit.jl/issues/3010
    * for some reason the ball model has has_hit_ground plot bugged, the step isnt shown in the plot idk why
- figure out how to end the simulation with an event?
- for the ball model, we cant pick a tspan we need to run it until the affect terminates the simulation (try SteadyStateProblem)
    * SteadyState was tried but failed with an instability RetCode idky