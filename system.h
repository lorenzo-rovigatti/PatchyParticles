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
void system_readColorsMax(char *namefile,int *max_species,int *max_colors);
void system_readColors(char *namefile,int **colorint,int *ncolorint,int **particle,int **color);
void system_readSpecies(char *nomefile,int n,int species,PatchyParticle *p,int *species_count);

#endif /* SYSTEM_H_ */
