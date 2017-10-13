# PatchyParticles

`PatchyParticles` is a code for perfoming Monte Carlo simulations of hard particles 
decorated with four patches modelled through the 
[Kern-Frenkel](http://www.sklogwiki.org/SklogWiki/index.php/Kern_and_Frenkel_patchy_model) 
pair interaction potential. The code is _educational_ in the sense that it is meant to 
be read as much as meant to be used.  

## Compilation

`PatchyParticles` does not require any dependencies but a C compiler, the standard C 
library and `make`. The default compiler, specified in the `makefile` file, is gcc. 
In order to compile the code it is sufficient to launch make within the main folder: 

	$ make 
	
Depending on the platform, it might be necessary to explicitly pass the makefile to
`make`. This can be done by issuing

	$ make -f makefile
	
The default behaviour is to compile `PatchyParticles` with full optimisations enabled
(that is, with the gcc options `-O3 -ffast-math -DNDEBUG`). The code can also be 
compiled with optimisations and the debug flags `-g2 -pg -ggdb3` (`make g=1`) or
with debug flags only (`make dbg=1`).  

At the end of the compilation stage, two executables will be generated: `PatchyParticles`
and `generator`.

## Usage

The `PatchyParticles` program takes one mandatory argument, the input file. The syntax 
of the file is quite simple, as it is just a list of "key = value" lines. The Examples 
folder contains a simple input file with all the mandatory options.

The `generator` binary can be used to generate initial configurations of N particles 
at density ρ, with ρ < 0.7. It requires two arguments (an integer for N and a decimal
number for ρ). The program prints the configuration in a new file named `generated.rrr`.

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
* The `utils.c, utils.h` pair contains commonly-used functions to work with vectors
and matrices.
* The `system.c, system.h` pair contains the functions used to initialise and cleanup
the main data structure (`System`).
* The `output.c, output.h` pair contains the functions used to initialise and cleanup
the data structure responsible for printing the output (`Output`), as well as the 
functions that actually print most of the output.
* The `parse_input.c, parse_input.h` pair contains the data structures, logic and 
functions used by the code to parse the input file passed to `PatchyParticles`  
