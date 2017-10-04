/*
 * system.h
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "defs.h"
#include "parse_input.h"

void init_patches(LR_system *syst, LR_IO *IO, FILE *patch_file);
void init_system(input_file *input, LR_system *syst, LR_IO *IO);
void clean_system(LR_system *system);

void check_cells(LR_system *syst, LR_IO *IO);

#endif /* SYSTEM_H_ */
