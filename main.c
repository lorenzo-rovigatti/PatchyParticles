/*
 * Patchy2A9B.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "defs.h"
#include "output.h"
#include "MC.h"
#include "parse_input.h"
#include "neighs.h"
#include "system.h"
#include "utils.h"

int stop = 0;

void gbl_terminate(int arg) {
	signal(arg, SIG_DFL);
	fprintf(stderr, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	stop = 1;
}

int main(int argc, char *argv[]) {
	Output output_files;
	output_files.log = stderr;

	if(argc == 1) output_exit(&output_files, "Usage is %s input\n", argv[0]);

	// here we handle a few SIG* signals;
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	/**
	 * Initialisation of the data structures
	 */
	input_file input;
	loadInputFile(&input, argv[1]);
	if(input.state == ERROR) exit(1);
	System syst;
	output_init(&input, &output_files);
	system_init(&input, &syst, &output_files);
	MC_init(&input, &syst, &output_files);

	int i;

	syst.energy = 0;
	for(i = 0; i < syst.N; i++)	{
		syst.energy += energy(&syst, syst.particles + i);
	}
	syst.energy *= 0.5;
	if(syst.overlap == 1) output_exit(&output_files, "Initial configuration contains an overlap, aborting\n");

	double E = (syst.N > 0) ? syst.energy / syst.N : 0;
	output_log_msg(&output_files, "Initial energy: %lf\n", E);

	// TODO: to be moved somewhere else
	/**
	 * Load from the input and configuration files the initial step, number of total steps and "save every" options
	 */
	llint steps, start_from, curr_step, print_every, save_every;
	getInputLLInt(&input, "Steps", &steps, 1);
	getInputLLInt(&input, "Print_every", &print_every, 1);
	getInputLLInt(&input, "Save_every", &save_every, 1);
	int restart_step_counter = 1;
	getInputInt(&input, "Restart_step_counter", &restart_step_counter, 0);
	if(restart_step_counter) {
		start_from = 0;
	}
	else {
		char name[256];
		getInputString(&input, "Initial_conditions_file", name, 1);
		FILE *conf = fopen(name, "r");
		fscanf(conf, "%lld %*d %*f %*f %*f\n", &start_from);
		fclose(conf);
		steps += start_from;
	}
	for(curr_step = start_from; curr_step < steps && !stop; curr_step++) {
		if((curr_step % print_every) == 0) {
			output_print(&output_files, &syst, curr_step);

			if(curr_step != 0) {
				check_energy(&syst, &output_files);
				check_cells(&syst, &output_files);

				if(syst.ensemble == SUS) output_sus(&output_files, &syst, curr_step);
			}
		}

		if(curr_step > 0 && (curr_step % save_every) == 0) {
			char name[1024];
			sprintf(name, "%s/conf_%lld.rrr", output_files.configuration_folder, curr_step);
			output_save(&output_files, &syst, curr_step, name);
			output_save(&output_files, &syst, curr_step, output_files.configuration_last);
		}

		syst.do_ensemble(&syst, &output_files);
	}

	output_save(&output_files, &syst, curr_step, output_files.configuration_last);
	E = (syst.N > 0) ? syst.energy / syst.N : 0;
	printf("%lld %lf %lf", curr_step, E, syst.energy);
	if(syst.ensemble != 0) printf(" %lf", syst.N / syst.V);
	printf("\n");

	MC_free(&syst);
	system_free(&syst);
	output_free(&output_files);

	cleanInputFile(&input);

	return 0;
}
