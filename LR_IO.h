/*
 * LR_IO.h
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#ifndef LR_IO_H_
#define LR_IO_H_

#include <stdarg.h>

#include "defs.h"
#include "parse_input.h"

void init_IO(input_file *input, LR_IO *IO);
void clean_IO(LR_IO *IO);

void print_sus(LR_IO *IO, LR_system *syst, llint step);
void print_e_sus(LR_IO *IO, LR_system *syst, llint step);
void _print_conf(LR_IO *IO, LR_system *syst, llint step, char *name);
void print_output(LR_IO *IO, LR_system *syst, llint step);
void print_conf(LR_IO *IO, LR_system *syst, llint step);
void log_msg(LR_IO *IO, char *format, ...);
void die(LR_IO *IO, char *format, ...);
void die_stderr(char *format, ...);

#endif /* LR_IO_H_ */
