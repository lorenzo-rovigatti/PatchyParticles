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

void system_init(input_file *input, System *syst, Output *IO);
void system_free(System *system);

#endif /* SYSTEM_H_ */
