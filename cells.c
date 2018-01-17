/*
 * cells.c
 *
 *  Created on: 16 Oct 2017
 *      Author: lorenzo
 */

#include "cells.h"

#include "defs.h"
#include "output.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>

void cells_init(System *syst, Output *output_files, double rcut) {
	syst->cells = (Cells *) malloc(sizeof(Cells));

	Cells *cells = syst->cells;
	int i;
	for(i = 0; i < 3; i++) {
		cells->N_side[i] = floor(syst->box[i] / rcut);
		if(cells->N_side[i] < 3) {
			cells->N_side[i] = 3;
			output_log_msg(output_files, "The size of the box along the %d-th dimension is too small, setting cells->N_side[%d] = 3\n", i, i);
		}
	}

	cells->N = cells->N_side[0] * cells->N_side[1] * cells->N_side[2];
	cells->heads = malloc(sizeof(PatchyParticle *) * cells->N);
	cells->next = malloc(sizeof(PatchyParticle *) * syst->N_max);
	output_log_msg(output_files, "Cells per side: (%d, %d, %d), total: %d\n", cells->N_side[0], cells->N_side[1], cells->N_side[2], cells->N);

	for(i = 0; i < cells->N; i++) cells->heads[i] = NULL;
	for(i = 0; i < syst->N_max; i++) cells->next[i] = NULL;
}

int cells_fill_and_get_idx_from_particle(System *syst, PatchyParticle *p, int idx[3]) {
	return cells_fill_and_get_idx_from_vector(syst, p->r, idx);
}

int cells_fill_and_get_idx_from_vector(System *syst, vector r, int idx[3]) {
	idx[0] = (int) ((r[0] / syst->box[0] - floor(r[0] / syst->box[0])) * (1. - DBL_EPSILON) * syst->cells->N_side[0]);
	idx[1] = (int) ((r[1] / syst->box[1] - floor(r[1] / syst->box[1])) * (1. - DBL_EPSILON) * syst->cells->N_side[1]);
	idx[2] = (int) ((r[2] / syst->box[2] - floor(r[2] / syst->box[2])) * (1. - DBL_EPSILON) * syst->cells->N_side[2]);

	return (idx[0] * syst->cells->N_side[1] + idx[1]) * syst->cells->N_side[2] + idx[2];
}

void cells_fill(System *syst) {
	int i, ind[3];
	Cells *cells = syst->cells;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;

		int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
		cells->next[p->index] = cells->heads[cell_index];
		cells->heads[cell_index] = p;
		p->cell = cell_index;
		p->cell_old = cell_index;
	}
}

void cells_check(System *syst, Output *output_files) {
	int i, counter;
	for(i = 0, counter = 0; i < syst->cells->N; i++) {
		PatchyParticle *p = syst->cells->heads[i];
		while(p != NULL) {
			p = syst->cells->next[p->index];
			counter++;
		}
	}

	if(counter != syst->N) output_exit(output_files, "\nThere are %d particles in cells, there should be %d\n", counter, syst->N);
}

void cells_free(Cells *cells) {
	free(cells->heads);
	free(cells->next);
	free(cells);
}
