/*
 * system.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "MC.h"
#include "LR_system.h"
#include "LR_IO.h"
#include "neighs.h"
#include "utils.h"

void init_cells(LR_system *syst, LR_IO *IO) {
	LR_cells *cells = &syst->cells;
	cells->N_side = floor(syst->L / (1. + syst->kf_delta_aa));
	if(cells->N_side < 3) {
		cells->N_side = 3;
		log_msg(IO, "Box side is too small, setting cells.N_side = 3\n");
	}
	cells->N = cells->N_side * cells->N_side * cells->N_side;
	cells->heads = malloc(sizeof(PatchyParticle *) * cells->N);
	log_msg(IO, "Cells per side: %d, total: %d\n", cells->N_side, cells->N);

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

void init_patches(LR_system *syst, LR_IO *IO, FILE *patch_file) {
	int i, j;

	if(patch_file == NULL) {
		switch(syst->n_patches) {
		case 1:
			set_vector(syst->base_patches[0], 0, 0, 1);
			break;
		case 2:
			set_vector(syst->base_patches[0], 0, 0, 1);
			set_vector(syst->base_patches[1], 0, 0, -1);
			break;
		case 3: {
			double cos30 = cos(M_PI / 6.);
			double sin30 = sin(M_PI / 6.);

			set_vector(syst->base_patches[0], 0, 1, 0);
			set_vector(syst->base_patches[1], cos30, -sin30, 0);
			set_vector(syst->base_patches[2], -cos30, -sin30, 0);
			break;
		}
		case 4: {
			double half_isqrt3 = 0.5 / sqrt(3);
			set_vector(syst->base_patches[0], -half_isqrt3, -half_isqrt3,  half_isqrt3);
			set_vector(syst->base_patches[1], half_isqrt3, -half_isqrt3, -half_isqrt3);
			set_vector(syst->base_patches[2], half_isqrt3,  half_isqrt3,  half_isqrt3);
			set_vector(syst->base_patches[3], -half_isqrt3,  half_isqrt3, -half_isqrt3);
			break;
		}
		default:
			die(IO, "Unsupported number of patches. Were you trying to use your own patch_file?\n");
		}
	}
	else {
		for(i = 0; i < syst->n_patches; i++) {
			// adding to j just to get rid of a warning
			j += fscanf(patch_file, "%lf %lf %lf\n", syst->base_patches[i], syst->base_patches[i] + 1, syst->base_patches[i] + 2);
		}
	}

	for(i = 0; i < syst->n_patches; i++) normalize(syst->base_patches[i]);

	// check that patches do not overlap with each other
	double theta_limit = 2*acos(syst->kf_cosmax_aa);
	for(i = 0; i < syst->n_patches; i++) {
		for(j = i+1; j < syst->n_patches; j++) {
			double theta = acos(SCALAR(syst->base_patches[i], syst->base_patches[j]));
			if(theta < theta_limit) die(IO, "Patches %d and %d are overlapping (the angle between them is %lf, while it should not be smaller than %lf)\n", i, j, theta, theta_limit);
		}
	}

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

void init_system(input_file *input, LR_system *syst, LR_IO *IO) {
	int res, i;

	getInputInt(input, "Dynamics", &syst->dynamics, 1);
	getInputInt(input, "Ensemble", &syst->ensemble, 1);

	getInputDouble(input, "DispDelta", &syst->disp_delta, 1);
	getInputDouble(input, "OrDelta", &syst->or_delta, 1);
	getInputDouble(input, "Temperature", &syst->T, 1);

	char name[256];
	getInputString(input, "Initial_conditions_file", name, 1);
	if(getInputInt(input, "Seed", &syst->seed, 0) == KEY_NOT_FOUND) {
		syst->seed = time(NULL);
		log_msg(IO, "Using seed %d\n", syst->seed);
	}
	srand48(syst->seed);

	FILE *conf = fopen(name, "r");
	if(conf == NULL) {
		int make_initial = 0;
		getInputInt(input, "make_initial", &make_initial, 0);
		if(!make_initial) {
			log_msg(IO, "The initial configuration can be generated putting 'make_initial = 1' in the input file\n");
			die(IO, "Initial_conditions_file '%s' is not readable\n", name);
		}

		// very slow way of making a configuration, not recommended
		LR_system new_syst;
		new_syst.N = 0;
		new_syst.ensemble = 0;
		new_syst.dynamics = 0;
		new_syst.kf_delta_aa = 0;
		new_syst.n_patches = 0;
		new_syst.base_patches = NULL;

		getInputDouble(input, "box_size", &new_syst.L, 1);
		init_cells(&new_syst, IO);
		getInputInt(input, "initial_N", &new_syst.N, 1);
		new_syst.particles = malloc(new_syst.N * sizeof(PatchyParticle));
		for(i = 0; i < new_syst.N; i++) {
			PatchyParticle *p = new_syst.particles + i;
			p->patches = NULL;
		}
		log_msg(IO, "Creating the configuration (N = %d, L = %lf), this could take some time... ", new_syst.N, new_syst.L);
		make_initial_conf(&new_syst, IO, name);
		conf = fopen(name, "r");
		log_msg(IO, "done\n");

		clean_system(&new_syst);
	}

	res = fscanf(conf, "%*d %*d %d %*d %*d\n", &syst->N);
	res += fscanf(conf, "%lf %*f %*f %*f %*f %*f\n", &syst->L);
	if(res != 2) die(IO, "The initial configuration file '%s' is empty or its headers are malformed\n", name);

	getInputDouble(input, "kf_delta_aa", &syst->kf_delta_aa, 1);
	getInputDouble(input, "kf_cosmax_aa", &syst->kf_cosmax_aa, 1);
	getInputInt(input, "n_patches", &syst->n_patches, 1);

	syst->kf_sqr_rcut_aa = SQR(1. + syst->kf_delta_aa);

	int is_iso = 0;
	getInputInt(input, "isotropic", &is_iso, 0);
	if(is_iso) {
		syst->iso_epsilon = 1.*syst->n_patches/12.;
		getInputDouble(input, "isotropic_range", &syst->iso_range, 1);
		syst->iso_sqr_range = SQR(syst->iso_range);
		if(syst->iso_sqr_range < syst->kf_sqr_rcut_aa) die(IO, "The range of the isotropic attraction (%lf) should be larger than the patchy one (%lf)", sqrt(syst->iso_sqr_range), sqrt(syst->kf_sqr_rcut_aa));
		log_msg(IO, "Isotropic interaction parameters: epsilon = %lf, range = %lf\n", syst->iso_epsilon, syst->iso_range);
	}
	else syst->iso_epsilon = 0.;

	log_msg(IO, "Patch parameters: n_patches = %d, cosmax = %lf, delta = %lf\n", syst->n_patches, syst->kf_cosmax_aa, syst->kf_delta_aa);

	syst->base_patches = malloc(sizeof(vector) * syst->n_patches);

	if(syst->ensemble != 0) {
		syst->substeps = syst->N*10;
		if(syst->substeps == 0) syst->substeps = 10;
		getInputDouble(input, "Activity", &syst->z, 1);
		switch(syst->ensemble) {
		case 1:
			getInputInt(input, "Number_colloids", &syst->N_max, 1);
			syst->N_min = 0;
			syst->N_max *= 5;
			break;
		case 6:
			getInputDouble(input, "Umbrella_sampling_energystep", &syst->SUS_e_step, 1);
			/* no break */
		case 3:
			getInputInt(input, "Umbrella_sampling_min", &syst->N_min, 1);
			getInputInt(input, "Umbrella_sampling_max", &syst->N_max, 1);
			if(syst->N < syst->N_min) die(IO, "Number of particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->N, syst->N_min);
			if(syst->N > syst->N_max) die(IO, "Number of particles %d is larger than Umbrella_sampling_max (%d)\n", syst->N, syst->N_max);
			break;
		default:
			die(IO, "Unsupported ensemble '%d'\n", syst->ensemble);
			break;
		}

		if(syst->ensemble == 3) {
			syst->SUS_hist = calloc(syst->N_max - syst->N_min + 1, sizeof(llint));
		}
		if(syst->ensemble == 6) {
			double max_E = syst->n_patches / 2.;
			syst->SUS_e_bins = (int) (max_E * syst->N_max / syst->SUS_e_step);
			syst->SUS_e_step = max_E * syst->N_max / syst->SUS_e_bins;

			int n_dim = syst->N_max - syst->N_min + 1;
			syst->SUS_e_hist = malloc(n_dim*sizeof(llint *));
			for(i = 0; i < n_dim; i++) syst->SUS_e_hist[i] = calloc(syst->SUS_e_bins, sizeof(llint));
		}
	}
	else {
		syst->N_max = syst->N;
		syst->substeps = syst->N;
	}

	syst->V = syst->L * syst->L * syst->L;
	syst->particles = malloc(syst->N_max * sizeof(PatchyParticle));
	syst->energy = 0;
	syst->overlap = 0;

	// this should be initialised even though we carry out regular non-avb simulations
	// since single_cluster simulations need it
	// pi * (delta^3 + 3*delta^2 + 3*delta) * (1 - cosmax)^2 / 3
	syst->avb_vin = syst->n_patches*syst->n_patches * (M_PI*(syst->kf_delta_aa*syst->kf_delta_aa*syst->kf_delta_aa + 3.*SQR(syst->kf_delta_aa) +
			3.*syst->kf_delta_aa) * SQR(1. - syst->kf_cosmax_aa)/3.);
	log_msg(IO, "Vavb = %lf\n", syst->avb_vin);
	if(syst->dynamics != 0) {
		syst->avb_vout = syst->V - syst->avb_vin;

		syst->avb_p = 0.5;
		getInputDouble(input, "avb_p", &syst->avb_p, 0);
	}

	syst->use_avb_in_in = 0;
	getInputInt(input, "avb_in_in", &syst->use_avb_in_in, 0);
	if(syst->use_avb_in_in) log_msg(IO, "AVB in-in enabled\n");

	char patch_file_name[512];
	FILE *patch_file = NULL;
	if(getInputString(input, "patch_file", patch_file_name, 0) == KEY_FOUND) {
		patch_file = fopen(patch_file_name, "r");
		if(patch_file == NULL) die(IO, "patch_file '%s' is not readable", patch_file_name);
	}

	init_patches(syst, IO, patch_file);
	if(patch_file != NULL) fclose(patch_file);
	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		p->index = i;
		p->type = TYPE_NO;
		p->n_aa = 0;

		p->patches = malloc(sizeof(vector) * syst->n_patches);
	}

	i = 0;
	vector p1, p2, p3;
	char myline[512];
	char *s_res = fgets(myline, 512, conf);
	while(s_res != NULL && strcmp(myline, "BONDS\n")) {
		sscanf(myline, "%lf %lf %lf\n", p1, p1 + 1, p1 + 2);
		res = fscanf(conf, "%lf %lf %lf\n", p2, p2 + 1, p2 + 2);

		PatchyParticle *p = syst->particles + i;
		res = fscanf(conf, "%lf %lf %lf\n", p->r, p->r+1, p->r+2);

		// normalize the orientation matrix
		normalize(p1);
		normalize(p2);
		cross(p1, p2, p3);

		// construct the orientation matrix
		memcpy(p->orient[0], p1, 3 * sizeof(double));
		memcpy(p->orient[1], p2, 3 * sizeof(double));
		memcpy(p->orient[2], p3, 3 * sizeof(double));
		gram_schmidt(p->orient[0], p->orient[1], p->orient[2]);

		int j;
		for(j = 0; j < syst->n_patches; j++) MATRIX_VECTOR_MULTIPLICATION(p->orient, syst->base_patches[j], p->patches[j]);

		i++;

		s_res = fgets(myline, 512, conf);
	}
	if(i != syst->N) die(IO, "Number of particles found in configuration (%d) is different from the value found in the header (%d)\n", i, syst->N);

	reset_counters(syst);

	init_cells(syst, IO);

	fclose(conf);
}

void clean_system(LR_system *syst) {
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
	if(syst->ensemble == 6) {
		for(i = 0; i < (syst->N_max - syst->N_min + 1); i++) free(syst->SUS_e_hist[i]);
		free(syst->SUS_e_hist);
	}
}

void check_cells(LR_system *syst, LR_IO *IO) {
	int i, counter;
	for(i = 0, counter = 0; i < syst->cells.N; i++) {
		PatchyParticle *p = syst->cells.heads[i];
		while(p != NULL) {
			p = p->next;
			counter++;
		}
	}

	if(counter != syst->N) die(IO, "\nThere are %d particles in cells, should be %d\n", counter, syst->N);
}
