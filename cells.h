/*
 * cells.h
 *
 *  Created on: 16 Oct 2017
 *      Author: lorenzo
 */

#ifndef CELLS_H_
#define CELLS_H_

typedef struct PatchyParticle PatchyParticle;
typedef struct Output Output;
typedef struct System System;

typedef struct Cells {
	int N_side[3];
	int N;
	PatchyParticle **heads;
	PatchyParticle **next;
} Cells;

void cells_init(System *syst, Output *IO, double rcut);
int cells_fill_and_get_idx(System *syst, PatchyParticle *p, int idx[3]);
void cells_fill(System *syst);
void cells_check(System *syst, Output *output_files);
void cells_free(Cells *cells);

#endif /* CELLS_H_ */
