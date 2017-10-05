/*
 * Patchy2A9B.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "defs.h"
#include "LR_IO.h"
#include "LR_system.h"
#include "MC.h"
#include "parse_input.h"
#include "neighs.h"
#include "utils.h"

int stop = 0;

void gbl_terminate(int arg) {
	signal(arg, SIG_DFL);
	fprintf(stderr, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	stop = 1;
}

int main(int argc, char *argv[]) {
	LR_IO IO;
	IO.log = stderr;

	if(argc == 1) die(&IO, "Usage is %s input\n", argv[0]);

	input_file input;
	loadInputFile(&input, argv[1]);
	if(input.state == ERROR) exit(1);

	// here we handle a few SIG* signals;
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	LR_system syst;

	init_IO(&input, &IO);
	init_system(&input, &syst, &IO);
	init_MC(&input, &syst, &IO);

	cleanInputFile(&input);

	int i;

	syst.energy = 0;
	for(i = 0; i < syst.N; i++)	{
		syst.energy += energy(&syst, syst.particles + i);
	}
	syst.energy *= 0.5;
	if(syst.overlap == 1) die(&IO, "Initial configuration contains an overlap, aborting\n");

	double E = (syst.N > 0) ? syst.energy / syst.N : 0;
	log_msg(&IO, "Initial energy: %lf\n", E);

	llint steps, curr_step, print_every, save_every;
	getInputLLInt(&input, "steps", &steps, 1);
	getInputLLInt(&input, "print_every", &print_every, 1);
	getInputLLInt(&input, "save_every", &save_every, 1);
	for(curr_step = 0; curr_step < steps && !stop; curr_step++) {
		if((curr_step % print_every) == 0) {
			print_output(&IO, &syst, curr_step);

			if(curr_step != 0) {
				log_msg(&IO, "Checking coherence...");
				check_energy(&syst, &IO);
				check_cells(&syst, &IO);
				log_msg(&IO, "done\n");

				if(syst.ensemble == 3) print_sus(&IO, &syst, curr_step);
				if(syst.ensemble == 6) print_e_sus(&IO, &syst, curr_step);
			}
		}

		if(curr_step > 0 && (curr_step % save_every) == 0) {
			print_conf(&IO, &syst, curr_step);
		}

		MC_sweep(&syst, &IO);
	}

	_print_conf(&IO, &syst, curr_step, IO.configuration_last);
	E = (syst.N > 0) ? syst.energy / syst.N : 0;
	printf("%lld %lf %lf", curr_step, E, syst.energy);
	if(syst.ensemble != 0) printf(" %lf", syst.N / syst.V);
	printf("\n");

	clean_system(&syst);
	clean_IO(&IO);

	return 0;
}
