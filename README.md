# MC_carbon_miswired_magnet
Simulation package for miswired magnet in Hall A 
This is modify version of SAMC from A. Deur.
This simulation is for the mis-wired magnet for 1st period experiment E97-110 at Hall A JLab.
Changes:
1. 4 forward matrices are used to transport particle from target to focal plane.
These matrices are obtained from fitting experimental data sieve slit holes to the survey information.
2. There is no apperture or any magnet inside the code.
Particle from target to focal plane by the transport matrices. There is no lost all the way.
Only loss from radiative correction if it is turn ON.

Input files:
Coefficients for the forward matrices: x, y, theta, phi. (x_coef,..ect)

The main code is mce97110h_c12.f

Input file c12.inp: contain information about:
Number of trial, beam energy, momentum setting(p0),angle, momentum range (+-x = 2x), theta range, phi range, beam x position (raster), beam y position (raster), target length, xxx (target orientation?), phase space option, total incoming radiation length (RL), outgoing radiation length (RL), average incoming RL density, average outgoing RL density.

Input cut file cutc12.inp: this file is not useful. But I leave it there in case need it later.

run.sh: this file to run the simulation many times to get enough statistic.
