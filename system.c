/*
 * system.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "system.h"

#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "MC.h"
#include "output.h"
#include "neighs.h"
#include "utils.h"

void _init_cells(System *syst, Output *IO) {
	Cells *cells = &syst->cells;
	cells->N_side = floor(syst->L / (1. + syst->kf_delta));
	if(cells->N_side < 3) {
		cells->N_side = 3;
		output_log_msg(IO, "Box side is too small, setting cells.N_side = 3\n");
	}
	cells->N = cells->N_side * cells->N_side * cells->N_side;
	cells->heads = malloc(sizeof(PatchyParticle *) * cells->N);
	output_log_msg(IO, "Cells per side: %d, total: %d\n", cells->N_side, cells->N);

	int i, ind[3];
	for(i = 0; i < cells->N; i++) cells->heads[i] = NULL;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;

		ind[0] = (int) ((p->r[0] / syst->L - floor(p->r[0] / syst->L)) * (1. - DBL_EPSILON) * cells->N_side);
		ind[1] = (int) ((p->r[1] / syst->L - floor(p->r[1] / syst->L)) * (1. - DBL_EPSILON) * cells->N_side);
		ind[2] = (int) ((p->r[2] / syst->L - floor(p->r[2] / syst->L)) * (1. - DBL_EPSILON) * cells->N_side);

		int cell_index = (ind[0] * cells->N_side + ind[1]) * cells->N_side + ind[2];
		p->next = cells->heads[cell_index];
		cells->heads[cell_index] = p;
		p->cell = cell_index;
		p->cell_old = cell_index;
	}
}

void _init_tetrahedral_patches(System *syst, Output *IO) {
	syst->n_patches = 4;
	syst->base_patches = malloc(sizeof(vector) * syst->n_patches);
	double half_isqrt3 = 0.5 / sqrt(3);
	set_vector(syst->base_patches[0], -half_isqrt3, -half_isqrt3,  half_isqrt3);
	set_vector(syst->base_patches[1], half_isqrt3, -half_isqrt3, -half_isqrt3);
	set_vector(syst->base_patches[2], half_isqrt3,  half_isqrt3,  half_isqrt3);
	set_vector(syst->base_patches[3], -half_isqrt3,  half_isqrt3, -half_isqrt3);

	int i, j;
	for(i = 0; i < syst->n_patches; i++) normalize(syst->base_patches[i]);

	// now we need to initialize syst->base_orient
	// first we initialize my_orient as the identity matrix
	// and we get -syst->base_patches[0], because the
	// set_orientation_around_vector invert its first argument
	matrix my_orient;
	vector my_first_patch;
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			memset(my_orient[i], 0, 3 * sizeof(double));
			my_orient[i][i] = 1.;
		}
		my_first_patch[i] = -syst->base_patches[0][i];
	}
	// then we calculate the rotation matrix required to transform
	// the 0, 0, 1 vector to the syst->base_patches[0] one
	set_orientation_around_vector(my_first_patch, my_orient, 0);
	// and then we transpose that matrix to obtain the rotation
	// needed to transform the syst->base_patches[0] vector
	// into the 0, 0, 1 one
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) syst->base_orient[i][j] = my_orient[j][i];
	}
}

void system_init(input_file *input, System *syst, Output *IO) {
	int res, i;

	getInputInt(input, "Dynamics", &syst->dynamics, 1);
	getInputInt(input, "Ensemble", &syst->ensemble, 1);

	getInputDouble(input, "Disp_max", &syst->disp_max, 1);
	getInputDouble(input, "Theta_max", &syst->theta_max, 1);
	getInputDouble(input, "Temperature", &syst->T, 1);

	char name[256];
	getInputString(input, "Initial_conditions_file", name, 1);
	if(getInputInt(input, "Seed", &syst->seed, 0) == KEY_NOT_FOUND) {
		syst->seed = time(NULL);
		output_log_msg(IO, "Using seed %d\n", syst->seed);
	}
	srand48(syst->seed);

	FILE *conf = fopen(name, "r");
	if(conf == NULL) output_exit(IO, "Initial_conditions_file '%s' is not readable\n", name);

	res = fscanf(conf, "%*d %d %lf %*f %*f\n", &syst->N, &syst->L);
	if(res != 2) output_exit(IO, "The initial configuration file '%s' is empty or its headers are malformed\n", name);

	getInputDouble(input, "KF_delta", &syst->kf_delta, 1);
	getInputDouble(input, "KF_cosmax", &syst->kf_cosmax, 1);

	syst->kf_sqr_rcut = SQR(1. + syst->kf_delta);

	output_log_msg(IO, "Patch parameters: cosmax = %lf, delta = %lf\n", syst->kf_cosmax, syst->kf_delta);

	if(syst->ensemble != NVT) {
		getInputDouble(input, "Activity", &syst->z, 1);
		switch(syst->ensemble) {
		case GC:
			getInputInt(input, "GC_N_max", &syst->N_max, 1);
			syst->N_min = 0;
			break;
		case SUS:
			getInputInt(input, "Umbrella_sampling_min", &syst->N_min, 1);
			getInputInt(input, "Umbrella_sampling_max", &syst->N_max, 1);
			if(syst->N < syst->N_min) output_exit(IO, "Number of particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->N, syst->N_min);
			if(syst->N > syst->N_max) output_exit(IO, "Number of particles %d is larger than Umbrella_sampling_max (%d)\n", syst->N, syst->N_max);
			syst->SUS_hist = calloc(syst->N_max - syst->N_min + 1, sizeof(llint));
			break;
		default:
			output_exit(IO, "Unsupported ensemble '%d'\n", syst->ensemble);
			break;
		}
	}
	else {
		syst->N_max = syst->N;
	}

	syst->V = syst->L * syst->L * syst->L;
	syst->particles = malloc(syst->N_max * sizeof(PatchyParticle));
	syst->energy = 0;
	syst->overlap = 0;

	if(syst->dynamics != RTMC) {
		syst->avb_vin = syst->n_patches*syst->n_patches * (M_PI*(syst->kf_delta*syst->kf_delta*syst->kf_delta + 3.*SQR(syst->kf_delta) +
				3.*syst->kf_delta) * SQR(1. - syst->kf_cosmax)/3.);
		output_log_msg(IO, "Vavb = %lf\n", syst->avb_vin);
		syst->avb_vout = syst->V - syst->avb_vin;

		syst->avb_p = 0.5;
		getInputDouble(input, "avb_p", &syst->avb_p, 0);
	}

	_init_tetrahedral_patches(syst, IO);
	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		p->index = i;

		p->n_patches = syst->n_patches;
		p->patches = malloc(sizeof(vector) * p->n_patches);
		p->base_patches = syst->base_patches;
	}

	i = 0;
	vector p1, p2, p3;
	char myline[512];
	char *s_res = fgets(myline, 512, conf);
	while(s_res != NULL) {
		sscanf(myline, "%lf %lf %lf\n", p1, p1 + 1, p1 + 2);
		res = fscanf(conf, "%lf %lf %lf\n", p2, p2 + 1, p2 + 2);

		PatchyParticle *p = syst->particles + i;
		res = fscanf(conf, "%lf %lf %lf\n", p->r, p->r+1, p->r+2);

		// normalize the orientation matrix
		normalize(p1);
		normalize(p2);
		cross(p1, p2, p3);

		// construct the orientation matrix
		memcpy(p->orientation[0], p1, 3 * sizeof(double));
		memcpy(p->orientation[1], p2, 3 * sizeof(double));
		memcpy(p->orientation[2], p3, 3 * sizeof(double));
		gram_schmidt(p->orientation[0], p->orientation[1], p->orientation[2]);

		int j;
		for(j = 0; j < p->n_patches; j++) {
			MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[j], p->patches[j]);
		}

		i++;

		s_res = fgets(myline, 512, conf);
	}
	fclose(conf);
	if(i != syst->N) output_exit(IO, "Number of particles found in configuration (%d) is different from the value found in the header (%d)\n", i, syst->N);

	utils_reset_acceptance_counters(syst);

	_init_cells(syst, IO);
}

void system_free(System *syst) {
	int i;
	if(syst->base_patches != NULL) free(syst->base_patches);
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		if(p->patches != NULL) {
			free(p->patches);
		}
	}
	free(syst->particles);
	free(syst->cells.heads);
	if(syst->ensemble == 3) free(syst->SUS_hist);
}

void check_cells(System *syst, Output *IO) {
	int i, counter;
	for(i = 0, counter = 0; i < syst->cells.N; i++) {
		PatchyParticle *p = syst->cells.heads[i];
		while(p != NULL) {
			p = p->next;
			counter++;
		}
	}

	if(counter != syst->N) output_exit(IO, "\nThere are %d particles in cells, should be %d\n", counter, syst->N);
}
