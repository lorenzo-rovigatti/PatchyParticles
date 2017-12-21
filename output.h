/*
 * output.h
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#ifndef output_H_
#define output_H_

#include <stdarg.h>

#include "defs.h"
#include "parse_input.h"

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
