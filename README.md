# VanillaCollisionSimulation

This code implements a collision simulation using the Lennard-Jones potential.

The following example shows the collision of two solid bodies that are relatively soft. The color of each particle encodes its velocity. Using the velocity as a marker for the particles' color allows to see the shockwaves that propagate through the bodies. It is interesting to note that the shock waves propagate at a significantly higher velocity in the medium compared to the velocity of the projectile.

In order to calculate the interaction of individual particles within a certain radius more efficiently, a tree algorithm was used to speed up the search for neighboring particles. This tree algorithm is based on code snippets of the book *Numerical simulation in molecular dynamics*.

The following two examples used the parameters shown below. For the second simulation the particle's mass of the smaller body was 100 times larger compared to the particles of the larger body.

```cpp
N_1 = 10000
N_2 = 100 // 64 for the second simulation
EPSILON = 5
SIGMA = 1
R_CUT = 2.5 * SIGMA
dt = 0.00005
PARTICLE_DISTANCE = 1.122462 * SIGMA // 2^(1/6) * SIGMA
```

**Example 1**

<p align="center">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/result_1.gif" height="500">
</p>

<div align="center">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-0.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-1.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-2.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-3.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-4.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-5.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-6.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-7.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res-8.png" height="200">
</div>

**Example 2**

<p align="center">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/result_sim_2.gif" height="500">
</p>

<div align="center">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_1.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_2.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_3.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_4.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_5.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_6.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_7.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_8.png" height="200">
<img src="https://github.com/KaiFabi/VanillaCollisionSimulation/blob/master/results/res_sim_2_9.png" height="200">
</div>

Compile and run the program using

`gcc -O -Wall collision.c -o collision -lm` and `./collision`
