# Program Goal

 This program solves the coupled ODEs governing planetary motion for two masses around a fixed primary mass using a relatively standard 4th order Runge-Kutta algorithm. Included in the program are a few namelist files with extension '.namelist' that contain the two planetary masses, the primary masses, initial positions, initial velocities, and the total number of sampling points and length of time for the orbits. The purpose for including the namelist files is for the program to read these files and solve the system based on supplied initial conditions and parameters as well as to illustrate two situations: (1) When both planets orbit the primary independently and (2) when one planet orbits the other as a satellite or "moon". Note that the primary mass is assumed to be held fixed at the origin. In addition, the planets are confined to move in the x-y plane. Furthermore, the program also calculates both the total angular momentum as well as the total energy at each time step. Included with the program files is a python jupyter notebook with extension '.ipynb' which plots the trajectory of the orbits for the case of circular planetary motion and for the case of a planet-moon system. These jupyter notebook file also plots the total angular momentum and total energy as a function of time for both orbit trajectories as well. One notices that both quantities are conserved (at least to standard numerical accuracy).


# Directions for Usage

Ensure all files referred to in the contents section are present in your current directory. Using your terminal, navigate to the relevant directory and type "make" into your terminal. Press "Enter". This will compile all files and create an executable
named `planetary`. Type "./planetary file_name.namelist" into your terminal, with "file_name" being whichever namelist file you choose, and press enter. Note that there are two namelist files included. Proceed with typing "./planetary file_name.namelist" into your terminal and pressing enter to generate results for the second namelist file. The results will be printed to a file which is given a name as dictated in the corresponding namelist file. If the default values are used, the program will print results to a file named `default_results.dat`. If the `circular.namelist` file is used, the program will print the results to file named `circular_results.dat`. If the `planet_moon.namelist` file is used, the results will be printed to a file named `planetmoon_results.dat`. After these results files have been generated, open the jupyter notebook file `plots_analysis5.ipynb` and run the program to visualize the orbits and conserved quantities.


# Contents  

`read_write.f90` Reads the input file (or default values if desired), and writes the positions, velocities, energies, and angular momenta for each body to a file. The file name is dependent on the namelist file.   

`ode_solver.f90` Contains the 4th order Runge-Kutta subroutine, along with the interface allowing RK-4 to solve generalized functions. 

`mechanics.f90` Defines the coupled system of planetary motion ODES as well as calculates the total angular momenta and total energy at each time step.

`types.f90` Defines the argument types, integers, and reals used within the program.

`plots_analysis5.ipynb` Contains routines for plotting the results. Click run all to visualize planetary trajectories and conserved quantities over time.

`makefile` Provides program the instructions to compile all .f90 files and create the executable.    

`main.f90` Contains the main calls to run the program.  

`circular.namelist` Included namelist file containing parameters for two masses in circular orbit around a fixed primary mass.

`planet_moon.namelist` Included namelist file for a planet-moon system orbiting a fixed primary mass.