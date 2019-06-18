# VanillaCollisionSimulation

This code implements a collision simulation using the Lennard-Jones potential.

The following example shows the collision of two solid bodies that are relatively soft. The color of each particle encodes its velocity. Using the velocity as marker for the particles color allows to see the shockwave that propagates through the bodies. First along the surface and later towards the center of the solid body.

For the simulation the following parameters were used:

```cpp
N_1 = 10000
N_2 = 100
EPSILON = 5
SIGMA = 1
R_CUT = 2.5 * SIGMA
dt = 0.00005
PARTICLE_DISTANCE = (2)^(1/6) * SIGMA
```
