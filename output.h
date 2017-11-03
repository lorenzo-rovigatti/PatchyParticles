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
void output_free(Output *IO);

void output_sus(Output *IO, System *syst, llint step);
void output_bsus(Output *IO, System *syst, llint step);
void output_e_sus(Output *IO, System *syst, llint step);
void output_save(Output *IO, System *syst, llint step, char *name);
void output_print(Output *IO, System *syst, llint step);
void output_log_msg(Output *IO, char *format, ...);
void output_exit(Output *IO, char *format, ...);
void output_exit_stderr(char *format, ...);

#endif /* output_H_ */
