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
#include <assert.h>

#include "MC.h"
#include "output.h"
#include "utils.h"

void _init_base_orient(System *syst) {
	int i, j;
	// we initialize my_orient as the identity matrix
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
		for(j = 0; j < 3; j++) {
			syst->base_orient[i][j] = my_orient[j][i];
		}
	}
}

void _init_patches(input_file *input, System *syst, Output *output_files) {
	char patch_filename[256];
	FILE *patch_file;
	int i;

	double base_kf_delta, base_kf_cosmax;

	getInputDouble(input, "KF_delta", &base_kf_delta, 1);
	getInputDouble(input, "KF_cosmax", &base_kf_cosmax, 1);

	int use_patch_file = 0;
	getInputInt(input, "Use_patch_file", &use_patch_file, 0);

	if(use_patch_file) {
		if(getInputString(input, "Patch_file", patch_filename, 0) == KEY_NOT_FOUND) {
			output_exit(output_files, "No patch file specified\n");
		}

		patch_file = fopen(patch_filename, "r");
		if(patch_file == NULL) {
			output_exit(output_files, "Patch_file '%s' is not readable\n", patch_filename);
		}

		output_log_msg(output_files, "Initialising patches from '%s'\n", patch_filename);

		if(fscanf(patch_file, "%d\n", &syst->n_patches) != 1) {
			output_exit(output_files, "The first line of the patch file should contain the number of patches\n");
		}
	}
	else {
		syst->n_patches = 4;
	}

	syst->base_patches = malloc(sizeof(vector) * syst->n_patches);
	syst->kf_cosmax = malloc(sizeof(vector) * syst->n_patches);
	syst->kf_delta = malloc(sizeof(vector) * SQR(syst->n_patches));
	syst->kf_interaction_matrix = malloc(sizeof(double) * SQR(syst->n_patches));
	// by default all patch-patch interactions are allowed, their depth is 1 and their KF parameters are the base ones
	for(i = 0; i < SQR(syst->n_patches); i++) {
		syst->kf_interaction_matrix[i] = 1;
		syst->kf_delta[i] = base_kf_delta;
	}

	for(i = 0; i < syst->n_patches; i++) {
		syst->kf_cosmax[i] = base_kf_cosmax;
	}

	if(use_patch_file) {
		char myline[512];
		for(i = 0; i < syst->n_patches; i++) {
			char *s_res = fgets(myline, 512, patch_file);
			if(s_res == NULL) {
				output_exit(output_files, "The patch file should contain a line per patch, but the file seems to have only %d lines instead of %d\n", i + 1, syst->n_patches + 1);
			}
			if(sscanf(myline, "%lf %lf %lf\n", syst->base_patches[i], syst->base_patches[i] + 1, syst->base_patches[i] + 2) != 3) {
				output_exit(output_files, "Line n.%d does not contain three fields\n", i + 1);
			}
		}
		fclose(patch_file);

		// check for overrides to the interaction matrix and to the KF parameters
		int j;
		input_file patch_input_file;
		loadInputFile(&patch_input_file, patch_filename, 0);
		for(i = 0; i < syst->n_patches; i++) {
			double value;
			char key[512];

			for(j = 0; j < syst->n_patches; j++) {
				// patch-patch well depth
				sprintf(key, "epsilon[%d][%d]", i, j);
				if(getInputDouble(&patch_input_file, key, &value, 0) == KEY_FOUND) {
					syst->kf_interaction_matrix[i + syst->n_patches * j] = syst->kf_interaction_matrix[j + syst->n_patches * i] = value;
				}
				// patch-patch range
				sprintf(key, "KF_delta[%d][%d]", i, j);
				if(getInputDouble(&patch_input_file, key, &value, 0) == KEY_FOUND) {
					syst->kf_delta[P_IDX(i, j)] = syst->kf_delta[P_IDX(j, i)] = value;
				}
			}

			// patch angular width
			sprintf(key, "KF_cosmax[%d]", i);
			if(getInputDouble(&patch_input_file, key, &value, 0) == KEY_FOUND) {
				syst->kf_cosmax[i] = value;
			}
		}
		cleanInputFile(&patch_input_file);
	}
	else {
		double half_isqrt3 = 0.5 / sqrt(3);
		set_vector(syst->base_patches[0], -half_isqrt3, -half_isqrt3, half_isqrt3);
		set_vector(syst->base_patches[1], half_isqrt3, -half_isqrt3, -half_isqrt3);
		set_vector(syst->base_patches[2], half_isqrt3, half_isqrt3, half_isqrt3);
		set_vector(syst->base_patches[3], -half_isqrt3, half_isqrt3, -half_isqrt3);

		output_log_msg(output_files, "Patch parameters: cosmax = %lf, delta = %lf\n", syst->kf_cosmax[0], syst->kf_delta[0]);
	}

	int normalise_patches = 1;
	getInputInt(input, "Normalise_patches", &normalise_patches, 0);
	if(normalise_patches) {
		for(i = 0; i < syst->n_patches; i++) {
			normalize(syst->base_patches[i]);
		}
	}

	_init_base_orient(syst);
}

void system_init(input_file *input, System *syst, Output *output_files) {
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
		output_log_msg(output_files, "Using seed %d\n", syst->seed);
	}
	srand48(syst->seed);

	FILE *conf = fopen(name, "r");
	if(conf == NULL) {
		output_exit(output_files, "Initial_conditions_file '%s' is not readable\n", name);
	}

	res = fscanf(conf, "%*d %d %lf %lf %lf\n", &syst->N, syst->box, syst->box + 1, syst->box + 2);
	if(res != 4) {
		output_exit(output_files, "The initial configuration file '%s' is empty or its headers are malformed\n", name);
	}

	if(syst->ensemble == SUS || syst->ensemble == GC || syst->ensemble == BSUS) {
		getInputDouble(input, "Activity", &syst->z, 1);
		switch(syst->ensemble) {
		case GC:
			getInputInt(input, "GC_N_max", &syst->N_max, 1);
			syst->N_min = 0;
			break;
		case SUS:
			getInputInt(input, "Umbrella_sampling_min", &syst->N_min, 1);
			getInputInt(input, "Umbrella_sampling_max", &syst->N_max, 1);
			if(syst->N < syst->N_min) {
				output_exit(output_files, "Number of particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->N, syst->N_min);
			}

			if(syst->N > syst->N_max) {
				output_exit(output_files, "Number of particles %d is larger than Umbrella_sampling_max (%d)\n", syst->N, syst->N_max);
			}

			syst->SUS_hist = calloc(syst->N_max - syst->N_min + 1, sizeof(llint));
			break;
		case BSUS:
			getInputInt(input, "Umbrella_sampling_min", &syst->N_min, 1);
			getInputInt(input, "Umbrella_sampling_max", &syst->N_max, 1);
			if(syst->N < syst->N_min) {
				output_exit(output_files, "Number of particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->N, syst->N_min);
			}

			if(syst->N > syst->N_max) {
				output_exit(output_files, "Number of particles %d is larger than Umbrella_sampling_max (%d)\n", syst->N, syst->N_max);
			}

			int transition_size = 3 * (syst->N_max - syst->N_min + 1);
			int histogram_size = syst->N_max - syst->N_min + 1;

			syst->bsus_collect = calloc(transition_size, sizeof(double));
			syst->bsus_tm = calloc(transition_size, sizeof(double));
			syst->bsus_normvec = calloc(histogram_size, sizeof(double));
			syst->bsus_pm = calloc(histogram_size, sizeof(double));

			char bsus_name[256];
			int bsus_value = getInputString(input, "Initial_bsus_file", bsus_name, 0);
			if(bsus_value != KEY_NOT_FOUND) {
				FILE *bsus_file = fopen(bsus_name, "rb");

				if(bsus_file) {
					output_log_msg(output_files, "Reading initial BSUS collection matrix\n");

					char myline[512];
					int p = 0;
					char *s_res = fgets(myline, 512, bsus_file);
					while(s_res != NULL) {
						sscanf(myline, "%lf %lf %lf\n", syst->bsus_collect + 3 * p, syst->bsus_collect + 3 * p + 1, syst->bsus_collect + 3 * p + 2);
						output_log_msg(output_files, "%lf %lf %lf\n", syst->bsus_collect[3 * p], syst->bsus_collect[3 * p + 1], syst->bsus_collect[3 * p + 2]);

						p++;

						s_res = fgets(myline, 512, bsus_file);
					}

					assert(p == histogram_size);

					fclose(bsus_file);

					bsus_update_histo(syst);
				}
				else {
					output_log_msg(output_files, "No initial BSUS collection matrix found\n");
				}
			}
			else {
				output_log_msg(output_files, "No initial BSUS collection matrix declared\n");
			}

			break;

		default:
			output_exit(output_files, "Unsupported ensemble '%d'\n", syst->ensemble);
			break;
		}
	}
	else {
		syst->N_max = syst->N;
	}

	syst->Lx_move = 0;
	if(syst->ensemble == NPT) {
		getInputDouble(input, "rescale_factor_max", &syst->rescale_factor_max, 1);
		getInputDouble(input, "P", &syst->P, 1);
	}
	else {
		getInputInt(input, "Lx_move", &syst->Lx_move, 0);
		if(syst->Lx_move) {
			if(syst->box[1] != syst->box[2]) {
				output_exit(output_files, "Lx_move = 1 requires Ly = Lz\n");
			}

			getInputDouble(input, "Lx_change_max", &syst->Lx_change_max, 1);
			getInputDouble(input, "Lyz_min", &syst->Lyz_min, 1);
			getInputDouble(input, "Lyz_max", &syst->Lyz_max, 1);

			if(syst->Lyz_min >= syst->box[1]) {
				output_exit(output_files, "Lyz_min should be smaller than the box size\n");
			}

			if(syst->Lyz_max <= syst->box[1]) {
				output_exit(output_files, "Lyz_max should be larger than the box size\n");
			}

			output_log_msg(output_files, "Ly and Lz are allowed to vary between %lf and %lf\n", syst->Lyz_min, syst->Lyz_max);
		}
	}

	syst->V = syst->box[0] * syst->box[1] * syst->box[2];
	syst->particles = malloc(syst->N_max * sizeof(PatchyParticle));
	syst->energy = 0;
	syst->overlap = 0;

	syst->shoulder_height = syst->shoulder_width = 0.;
	getInputDouble(input, "Repulsive_shoulder_width", &syst->shoulder_width, 0);
	getInputDouble(input, "Repulsive_shoulder_height", &syst->shoulder_height, 0);
	syst->shoulder_rcut_sqr = SQR(1. + syst->shoulder_width);

	syst->shoulder_lr_model = 0;
	getInputInt(input, "Repulsive_LR_model", &syst->shoulder_lr_model, 0);

	_init_patches(input, syst, output_files);
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
		res = fscanf(conf, "%lf %lf %lf\n", p->r, p->r + 1, p->r + 2);

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
	if(i != syst->N) {
		output_exit(output_files, "Number of particles found in configuration (%d) is different from the value found in the header (%d)\n", i, syst->N);
	}

	utils_reset_acceptance_counters(syst);

	syst->r_cut = 1. + syst->shoulder_width;
	for(i = 0; i < SQR(syst->n_patches); i++) {
		double new_r_cut = 1. + syst->kf_delta[i];
		if(new_r_cut > syst->r_cut) {
			syst->r_cut = new_r_cut;
		}
	}
	syst->sqr_rcut = SQR(syst->r_cut);

	cells_init(syst, output_files, syst->r_cut);
	cells_fill(syst);
}

void system_free(System *syst) {
	cells_free(syst->cells);

	int i;
	if(syst->base_patches != NULL) {
		free(syst->base_patches);
		free(syst->kf_cosmax);
		free(syst->kf_delta);
		free(syst->kf_interaction_matrix);
	}
	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		if(p->patches != NULL) {
			free(p->patches);
		}
	}

	free(syst->particles);

	if(syst->ensemble == SUS) {
		free(syst->SUS_hist);
	}
	else if(syst->ensemble == BSUS) {
		free(syst->bsus_collect);
		free(syst->bsus_tm);
		free(syst->bsus_normvec);
		free(syst->bsus_pm);
	}
}
