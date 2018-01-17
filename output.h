/*
 * output.h
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#ifndef output_H_
#define output_H_

#include <stdio.h>
#include <stdarg.h>

typedef struct input_file input_file;
typedef struct System System;
typedef long long int llint;

typedef struct Output {
	llint start_from;
	llint save_every;
	llint print_every;
	int print_bonds;
	int restart_step_counter;
	char configuration_folder[512];
	char configuration_last[512];
	char sus_folder[512];
	char bonds_folder[512];
	FILE *log;
	FILE *energy;
	FILE *density;
	FILE *acc;
} Output;

void output_init(input_file *input, Output *output_files);
void output_free(Output *output_files);

void output_sus(Output *output_files, System *syst, llint step);
void output_save(Output *output_files, System *syst, llint step, char *name);
void output_print(Output *output_files, System *syst, llint step);
void output_print_bonds(Output *output_files, System *syst, char *name);
void output_log_msg(Output *output_files, char *format, ...);
void output_exit(Output *output_files, char *format, ...);
void output_exit_stderr(char *format, ...);

#endif /* output_H_ */
