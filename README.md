# PatchyParticles

`PatchyParticles` is a code for perfoming Monte Carlo simulations of hard particles 
decorated with four patches modelled through the 
[Kern-Frenkel](http://www.sklogwiki.org/SklogWiki/index.php/Kern_and_Frenkel_patchy_model) 
pair interaction potential. The code is _educational_ in the sense that it is meant to 
be read as much as meant to be used.  

## Compilation

`PatchyParticles` does not require any dependencies but a C compiler, the standard C 
library and GNU's `make`. The default compiler, specified in the `makefile` file, is 
gcc, but the code can be compiled with icc as well. The compilation has been tested 
with gcc >= 4.8 and icc >= 13.
In order to compile the code it is sufficient to launch make within the main folder: 

	$ make 

Depending on the platform, it might be necessary to explicitly pass the makefile to
`make`. This can be done by issuing

	$ make -f makefile

The code can be compiled with icc with the command 

	$ make -f makefile.icc
	
The default behaviour is to compile `PatchyParticles` with full optimisations enabled
(that is, with the options `-O3 -ffast-math -DNDEBUG` or `-O3 -fast -DNDEBUG` if the 
code is compiled with gcc or icc, respectively)). The code can also be compiled with 
optimisations and the debug flags `-g2 -pg -ggdb3` (`make g=1`) or with debug flags only 
(`make dbg=1`).

At the end of the compilation stage, two executables will be generated: `PatchyParticles`
and `generator`.

## Usage

The `PatchyParticles` program takes one mandatory argument, the input file. The syntax 
of the file is quite simple, as it is just a list of "key = value" lines. The Examples 
folder contains some simple input files with all the mandatory options.

The `generator` binary can be used to generate initial configurations of N particles 
at density ρ, with ρ < 0.7. It requires two arguments (an integer for N and a decimal
number for ρ). The program prints the configuration in a new file named `generated.rrr`.

## Input file

Here is a list of mandatory options. Please refer to the input files in the `Examples` 
folder for common values. 

### Simulation options

* `Dynamics = <int>`: use 0 for rototranslations, 1 for AVB and 2 for VMMC.
* `Ensemble = <int>`: use 0 for NVT, 1 for Grand Canonical (muVT), 3 for Successive 
Umbrella Sampling.
* `Temperature = <float>`: temperature of the simulation, in units of the patch-patch 
bond.
* `Steps = <int>`: length of the simulation, in Monte Carlo steps.
* `GC_N_max = <int>`: maximum number of particles, above which the simulation will be 
stopped. Meaningful only for Grand Canonical (or SUS) simulations

### Input/output options

* `Initial_conditions_file`: initial configuration.
* `Print_every = <int>`: output frequency for energy, density and acceptances.
* `Save_every = <int>`: output frequency for configurations.
* `Energy_file = <string>`: name of the output file for the energy. 
* `Density_file = <string>`: name of the output file for the density.
* `Configuration_last = <string>`: name of the file which stores the last configuration 
(printed with frequency `Save_every` and at the end of the simulation)
* `Confguration_folder`: name of the folder which will contain the configuration files

### Kern-Frenkel options

* `KF_delta = <float>`: Radial width of the Kern-Frenkel patches.
* `KF_cosmax = <float>`: Angular width of the Kern-Frenkel patches.

### Monte Carlo moves options

* `Disp_max = <float>`: maximum trial displacement for translations.
* `Theta_max = <float>`: maximum trial angular displacement for rotations, in radians.
* `vmmc_max_move = <float>`: maximum allowed displacement for VMMC moves: if a vmmc move 
attempts to move a particle for more than this value, the move will be rejected.
* `vmmc_max_cluster = <int>`: maximum cluster size for VMMC moves: if a vmmc move attempts 
to move more than this number of particles, the move will be rejected.

## Code organisation

* The `main.c` file contains calls to the initialisation functions, the calculation
of the initial energy, the main loop and the calls to the cleanup functions.
* The general data structures used throughout the code, as well as some useful macros,
are stored in the `defs.h` file.
* The `MC.c, MC.h` pair contains the main logic of the code. It manages the calculation
of the energy, the particle rototranslations, the simulation of the different ensembles,
_etc._ 
* The `avb.c, avb.h` and `vmmc.c, vmmc.h` pairs contain the data structures and 
functions pertaining to the Aggregation-Volume-Bias and Virtual-Move-Monte-Carlo moves
that can be optionally enabled to speed up the simulation efficiencies.
* The `cells.c, cells.h` pair contain the data structures and functions pertaining to the
linked-lists used to keep track of the list of neighbours of each particle.
* The `utils.c, utils.h` pair contains commonly-used functions to work with vectors
and matrices.
* The `system.c, system.h` pair contains the functions used to initialise and cleanup
the main data structure (`System`).
* The `output.c, output.h` pair contains the functions used to initialise and cleanup
the data structure responsible for printing the output (`Output`), as well as the 
functions that actually print most of the output.
* The `parse_input.c, parse_input.h` pair contains the data structures, logic and 
functions used by the code to parse the input file passed to `PatchyParticles`  
