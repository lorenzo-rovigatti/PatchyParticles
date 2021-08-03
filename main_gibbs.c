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
	Output output_files1,output_files2;
	output_files1.log = stderr;
	output_files2.log = stderr;

	if(argc != 3) output_exit(&output_files1, "%s [input box 1] [input box 2]\n", argv[0]);

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
	input_file input1,input2;
	loadInputFile(&input1, argv[1]);
	loadInputFile(&input2, argv[2]);
	if(input1.state == ERROR) exit(1);
	if(input2.state == ERROR) exit(1);

	System syst1,syst2;
	output_init(&input1, &output_files1);
	output_init(&input2, &output_files2);
	system_init(&input1, &syst1, &output_files1);
	system_init(&input2, &syst2, &output_files2);
	MC_init(&input1, &syst1, &output_files1);
	MC_init(&input2, &syst2, &output_files2);

	/**
	 * Compute the initial energy and check whether there are overlaps in the initial configuration
	 */
	int i;
	syst1.energy = 0;
	for(i = 0; i < syst1.N; i++)	{
		syst1.energy += MC_energy(&syst1, syst1.particles + i);
		if(syst1.overlap == 1) output_exit(&output_files1, "Initial configuration contains an overlap, aborting\n");
	}
	syst1.energy *= 0.5;

	double E = (syst1.N > 0) ? syst1.energy / syst1.N : 0;
	output_log_msg(&output_files1, "Initial energy: %lf\n", E);

	syst2.energy = 0;
	for(i = 0; i < syst2.N; i++)	{
		syst2.energy += MC_energy(&syst2, syst2.particles + i);
		if(syst2.overlap == 1) output_exit(&output_files2, "Initial configuration contains an overlap, aborting\n");
	}
	syst2.energy *= 0.5;

	E = (syst2.N > 0) ? syst2.energy / syst2.N : 0;
	output_log_msg(&output_files2, "Initial energy: %lf\n", E);



	/**
	 * Get the number of steps to be run from the input file 1
	 */
	llint steps,curr_step;
	getInputLLInt(&input1, "Steps", &steps, 1);
	steps += output_files1.start_from;

	/**
	 * Main loop
	 */


	for(curr_step = output_files1.start_from; curr_step < steps && !stop; curr_step++) {
		/**
		 * Print the output (energy, density, acceptance, etc.) every "print_every" steps
		 */



		if((curr_step % output_files1.print_every) == 0) {
			output_print(&output_files1, &syst1, curr_step);
			output_print(&output_files2, &syst2, curr_step);

			if(curr_step != 0) {
				cells_check(&syst1, &output_files1);
				MC_check_energy(&syst1, &output_files1);

				cells_check(&syst2, &output_files2);
				MC_check_energy(&syst2, &output_files2);

			}
		}

		/**
		 * Print the configuration every "save_every" steps
		 */
		if(curr_step > 0 && (curr_step % output_files1.save_every) == 0) {

			char name[1024];
			sprintf(name, "%s/conf_%lld.rrr", output_files1.configuration_folder, curr_step);
			output_save(&output_files1, &syst1, curr_step, name);
			output_save(&output_files1, &syst1, curr_step, output_files1.configuration_last);

			if(output_files1.save_also_as_mgl) {
				sprintf(name, "%s/conf_%lld.mgl", output_files1.configuration_folder, curr_step);
				output_save_to_mgl(&output_files1, &syst1, name);
				output_save_to_mgl(&output_files1, &syst1, "last.mgl");

			}

			sprintf(name, "%s/conf_%lld.rrr", output_files2.configuration_folder, curr_step);
			output_save(&output_files2, &syst2, curr_step, name);
			output_save(&output_files2, &syst2, curr_step, output_files2.configuration_last);

			if(output_files2.save_also_as_mgl) {
				sprintf(name, "%s/conf_%lld.mgl", output_files2.configuration_folder, curr_step);
				output_save_to_mgl(&output_files2, &syst2, name);
				output_save_to_mgl(&output_files2, &syst2, "last.mgl");

			}


		}

		/**
		 * Perform a Monte Carlo sweep
		 */
		do_GIBBS(&syst1,&syst2,&output_files1,&output_files2,curr_step);
		//syst1.do_ensemble(&syst1, &output_files1);
		//syst2.do_ensemble(&syst2, &output_files2);

	}

	/**
	 * Print the last configuration and the last line of the output
	 */
	output_save(&output_files1, &syst1, curr_step, output_files1.configuration_last);
	if(output_files1.save_also_as_mgl) {
		output_save_to_mgl(&output_files1, &syst1, "last.mgl");
	}
	output_save(&output_files2, &syst2, curr_step, output_files2.configuration_last);
	if(output_files2.save_also_as_mgl) {
		output_save_to_mgl(&output_files2, &syst2, "last.mgl");
	}
	double E1 = (syst1.N > 0) ? syst1.energy / syst1.N : 0;
	double E2 = (syst2.N > 0) ? syst2.energy / syst2.N : 0;
	printf("%lld %lf %lf", curr_step, E1, E2);
	if(syst1.ensemble != 0) printf(" %lf", syst1.N / syst1.V);
	if(syst2.ensemble != 0) printf(" %lf", syst2.N / syst2.V);
	printf("\n");

	/**
	 * Cleanup
	 */
	MC_free(&syst1);
	system_free(&syst1);
	output_free(&output_files1);
	cleanInputFile(&input1);

	MC_free(&syst2);
	system_free(&syst2);
	output_free(&output_files2);
	cleanInputFile(&input2);


	return 0;
}
