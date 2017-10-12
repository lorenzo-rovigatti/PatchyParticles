/*
 * defs.h
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#ifndef DEFS_H_
#define DEFS_H_

#define SQR(x) ((x) * (x))
#define SCALAR(x, y) ((x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2])

#define MAX_E 2.5

#define N_MOVES 6
#define ROTO_TRASL 0
#define AVB 1
#define ADD 2
#define REMOVE 3
#define AVB_IN_IN 4
#define MOVE_VMMC 5

#define NVT 0
#define GC 1
#define SUS 3

#define RTMC 0
#define VMMC 1
#define AVBMC 2

#define TYPE_NO -1
#define TYPE_MON 0
#define TYPE_END 1

#define NO_BOND 0
#define ISO_BOND 1
#define PATCH_BOND 2

#define NO_NEIGH NULL

// Bad-looking but excellent performance-wise
#define MATRIX_VECTOR_MULTIPLICATION(m, v, result) {\
	(result)[0] = (m)[0][0]*(v)[0] + (m)[0][1]*(v)[1] + (m)[0][2]*(v)[2];\
	(result)[1] = (m)[1][0]*(v)[0] + (m)[1][1]*(v)[1] + (m)[1][2]*(v)[2];\
	(result)[2] = (m)[2][0]*(v)[0] + (m)[2][1]*(v)[1] + (m)[2][2]*(v)[2];\
}

#include <stdio.h>
#include <complex.h>
#include <assert.h>
#include <signal.h>

typedef double vector[3];
typedef double matrix[3][3];
typedef long long int llint;

typedef struct PatchyParticle {
	vector r, r_old;
	matrix orientation, orientation_old;
	int index;

	int n_patches;
	vector *patches, *base_patches;

	int cell, cell_old;
	struct PatchyParticle *next;
} PatchyParticle;

typedef struct Cells {
	int N_side;
	int N;
	PatchyParticle **heads;
	// TODO: implement this
	PatchyParticle **next;
} Cells;

typedef struct Output {
	llint start_from;
	llint save_every;
	llint print_every;
	int restart_step_counter;
	char configuration_folder[512];
	char configuration_last[512];
	char sus_folder[512];
	FILE *log;
	FILE *energy;
	FILE *density;
	FILE *acc;
} Output;

typedef struct System {
	int N, N_min, N_max;
	double L;
	double V;
	double T;
	double z;
	double energy;
	int n_patches;

	// this matrix allows base_patches[0] to be transformed
	// to 0, 0, 1
	matrix base_orient;
	vector *base_patches;

	int dynamics;
	int ensemble;
	void (*do_dynamics)(struct System *, Output *);
	void (*do_ensemble)(struct System *, Output *);

	llint *SUS_hist;
	llint **SUS_e_hist;
	double SUS_e_step;
	int SUS_e_bins;

	int overlap;

	int tries[N_MOVES];
	int accepted[N_MOVES];

	double disp_max;
	double theta_max;

	double kf_delta, kf_cosmax, kf_sqr_rcut;

	Cells cells;

	int seed;
	PatchyParticle *particles;
} System;

#endif /* DEFS_H_ */
