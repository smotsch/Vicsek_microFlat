Program *MicroVic_flat*
=======================

**Table of Contents**
- [Installation](#installation)
- [Program execution](#program-execution)
- [Output](#output)
- [Graphics](#graphics)
- [Parameters](#parameters)
- [The program](#the-program)

## Installation

Execute the `makefile`:
```bash
	> make
```
The compiler creates the executable file `./MicroVic_flat` and move it to the directory `bin`. By default, the makefile uses the *gfortran* compiler. If you want to use
*ifort*, open the `makefile` and use the instruction written at the bottom. 

## Program execution

To execute the program, run the binary file:
```bash
	> ./MicroVic_flat
```
The program is going to read the parameters of the model in the files:
* `PARAMETER_MicroVic_flat.txt`: for most of the parameters (number of particles, their speed...)
* `PARAMETER_init.txt`: for the initial condition (uniform, Gaussian...)

The parameters are written in external files because we do not need to recompile when we
change one parameter by doing so. The problem with this method is that we have to write
the value of the parameters at a given line. We cannot add a comment line in the file
`PARAMETER_MicroVic_flat.txt`, otherwise the numbering is ruined.

During the execution, the program  `MicroVic_flat` is going to compute the positions
(`X`) and the velocity (`theta`) of the particles at each time step. In the
terminal, the program gives some information about the parameters used for the
simulation. At the end, it displays the computation time.

## Output

The program `MicroVic_flat` can store two types of data. First, it can save the
trajectories and the velocities of each particle over time. At each time step,
it creates in the directory 'data' the files:
* `particleX_******`     : x coordinate of the particles
* `particleY_******`     : y coordinate of the particles
* `particleTheta_******` : velocity angle of the particles,

with `******` a counter of time step. The program writes only if the parameter
`isTrajectorySave` is `True` (line 28 in `PARAMETER_MicroVic_flat.txt`).

The other output possibilities are the  *macroscopic* quantities (density,
flux, global distribution of velocity). Still in the repertory `data`, the
program can write the following files at each time step:
* `rho1Dx_*`, `u1Dx_*`, `v1Dx_*` : density and mean velocity in the direction `x`
* `rho2D_*`, `u2D_*`, `v2D_*`    : density and mean velocity in 2D
* `densTheta_*`                  : distribution of velocity angle

It requires to have strictly positive value for the following parameters:
* `dx`     : mesh size of the grid in `x`
* `dxy`    : same thing but in 2D (same mesh size in `x` and `y`)
* `dtheta `: mesh size of the grid for the distribution of angle

If one parameter is zero, the program neither computes or writes the
corresponding macroscopic quantity. The method used to compute the density is
the PIC method (Particle-In-Cell) of order 2 by default.


## Graphics

To display the results of the computations with **Octave** or **Python**
```bash
	> Display_particlesFlat2D.m   # plot the particles
	> Display_macroFlat2D.m       # plot the macro. density/velocity
```
in the folder `visualization`.

## Parameters

* `PARAMETER_MicroVic_flat.txt`
 * `N`          : number of particles
 * `c`          : speed of the particles
 * `nu`         : intensity of the relaxation to the mean velocity
 * `d`          : intensity of the noise
 * `R`          : radius of interaction of the particle
 * `Time`       : total time of the simulation
 * `dt`         : time step
 * `Lx`         : size of the domain in x
 * `Ly`         : same thing in y
 * `boundCond`  : boundary condition
 * `choiceInit` : choice for the initial condition
 * `isGrid`     : use or not the trick of the grid
 * `isInitRand` : initialize or not the random
 * `isTrajectorySave` : save or note the trajectories
 * `dx`         : size of the meshgrid for the distributions in 1D
 * `dxy`        : same thing but in 2D
 * `dtheta`     : same thing for the distribution of angle
 * `jumpPrint`  : to save only at certain time step
* `PARAMETER_init.txt`
 * `initCondX`     : choice for the initial condition for X
 * `xMean`,`yMean`,`xStd`,`yStd` : parameters for the Gaussian
 * `proportion`    : proportion of particles at the left side (Riemann problem)
 * `initCondTheta` : choice for the initial condition for theta
 * `thetaMean`,`thetaStd` : parameters for the Gaussian
 * `thetaL`,`thetaR`,`temperature` : parameters for the Riemann problem


## The program

The main program is the file *main_MicroVic_flat.f90*.
It uses different modules (defined in separated files):

| File                              | Description   |
| ----------------------------------|:-------------:|
| `toolkit`                         | contains the usual functions (AngleVec, RandNorm...)
| `input_output_MicroVic_flat`      | for the input/output (Lecture, FilePrint...)
| `boundary_MicroVic_flat`          | the effect of the wall (Wall)
| `initial_condition_MicroVic_flat` | to initialize with the proper initial conditions (InitCond)
| `grid_interaction`                | to use the super grid (CellNumber, ListAdd...)
| `interaction_MicroVic_flat`       | to compute the average velocity around a particle
| `stat_MicroVic_flat`              | to compute and save the macroscopic quantities (Moment2D...)


The architecture of the program is the following:
```bash
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%    main_MicroVic_flat                                                      %%
  %%      -> declaration of variables                                           %%
  %%      -> lecture of parameters ("Lecture")                                  %%
  %%      -> initialization of variables ("InitCond")                           %%
  %%                                                                            %%
  %%      -> loop in time                                                       %%
  %%        1) computation of the average velocity Ω                            %%
  %%        2) new positions and velocities                                     %%
  %%        3) update of the grid (if necessary)                                %%
  %%        4) write trajectories and/or densities ("FilePrint")                %%
  %%      -<                                                                    %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
There is no Jedi mind tricks in the program. It supposes to be easy understandable.
There might be two points that needs some explanations:
* A structure `PARAM_MicroVic_flat` is used in order to save all the parameters of
the model in only one variable (just called `P`). This structure is defined in
the file `input_output_MicroVic_flat.f90`. This avoid to write 36 arguments each time a
subroutine is called.
* The grid method (also called the **Verlet list** by Wikipedia) consists to
allocate a  number at each particle depending on its position on a **virtual
grid** (`PosGrid`). Therefore, when we look for the neighbors of a particle,
we only have to   check the particles within a nearby square. The time-saving
is spectacular.


## Legal notice information

 This (modest) program is distributed under the GNU GPL license version 2. More
information are available in the file COPYING.txt. For any information or bugs,
please contact me at: `smotsch[at]asu.edu`.
