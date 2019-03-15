/*
 * Patchy2A9B.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "defs.h"
#include "MC.h"
#include "output.h"
#include "parse_input.h"
#include "system.h"
#include "utils.h"

int stop = 0;
void gbl_terminate(int arg) {
	/**
	 * The next time the signal is intercepted, the default handler will be invoked
	 * in order to avoid making the program unkillable.
	 */
	signal(arg, SIG_DFL);
	fprintf(stderr, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	stop = 1;
}

int main(int argc, char *argv[]) {
	Output output_files;
	output_files.log = stderr;

	if(argc == 1) output_exit(&output_files, "Usage is %s input\n", argv[0]);

	/**
	 * Here we handle a few SIG* signals: whenever we intercept one of these signals
	 * the program will quit the main loop and die as gracefully as possible.
	 */
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	/**
	 * Initialise the data structures
	 */
	input_file input;
	loadInputFile(&input, argv[1]);
	if(input.state == ERROR) exit(1);
	System syst;
	output_init(&input, &output_files);
	system_init(&input, &syst, &output_files);
	MC_init(&input, &syst, &output_files);

	/**
	 * Compute the initial energy and check whether there are overlaps in the initial configuration
	 */
	int i;
	syst.energy = 0;
	for(i = 0; i < syst.N; i++)	{
		syst.energy += MC_energy(&syst, syst.particles + i);
		if(syst.overlap == 1) output_exit(&output_files, "Initial configuration contains an overlap, aborting\n");
	}
	syst.energy *= 0.5;

	double E = (syst.N > 0) ? syst.energy / syst.N : 0;
	output_log_msg(&output_files, "Initial energy: %lf\n", E);

	/**
	 * Get the number of steps to be run from the input file
	 */
	llint steps, curr_step;
	getInputLLInt(&input, "Steps", &steps, 1);
	steps += output_files.start_from;

	/**
	 * Main loop
	 */


	for(curr_step = output_files.start_from; curr_step < steps && !stop; curr_step++) {
		/**
		 * Print the output (energy, density, acceptance, etc.) every "print_every" steps
		 */
		if((curr_step % output_files.print_every) == 0) {
			output_print(&output_files, &syst, curr_step);

			if(curr_step != 0) {
				MC_check_energy(&syst, &output_files);
				cells_check(&syst, &output_files);

				//if(syst.ensemble == SUS) output_sus(&output_files, &syst, curr_step);

				if(syst.ensemble == BSUS) bsus_update_histo(&syst);
			}
		}

		/**
		 * Print the configuration every "save_every" steps
		 */
		if(curr_step > 0 && (curr_step % output_files.save_every) == 0) {

			char name[1024];
			sprintf(name, "%s/conf_%lld.rrr", output_files.configuration_folder, curr_step);
			output_save(&output_files, &syst, curr_step, name);
			output_save(&output_files, &syst, curr_step, output_files.configuration_last);

			if(output_files.save_also_as_mgl) {
				sprintf(name, "%s/conf_%lld.mgl", output_files.configuration_folder, curr_step);
				output_save_to_mgl(&output_files, &syst, name);
				output_save_to_mgl(&output_files, &syst, "last.mgl");

			}
			if(syst.ensemble == BSUS)
			{
				output_bsus(&output_files, &syst, curr_step);
			}
		}

		/**
		 * Perform a Monte Carlo sweep
		 */
		syst.do_ensemble(&syst, &output_files);
	}

	/**
	 * Print the last configuration and the last line of the output
	 */
	output_save(&output_files, &syst, curr_step, output_files.configuration_last);
	if(output_files.save_also_as_mgl) {
		output_save_to_mgl(&output_files, &syst, "last.mgl");
	}
	E = (syst.N > 0) ? syst.energy / syst.N : 0;
	printf("%lld %lf %lf", curr_step, E, syst.energy);
	if(syst.ensemble != 0) printf(" %lf", syst.N / syst.V);
	printf("\n");

	/**
	 * Cleanup
	 */
	MC_free(&syst);
	system_free(&syst);
	output_free(&output_files);
	cleanInputFile(&input);

	return 0;
}
