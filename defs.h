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

#define RESOLVE_OVERLAP_LIMIT 1000000
#define JOIN_CLUSTERS_LIMIT 100000

#define MAX_E 2.5

#define N_MOVES 6
#define ROTO_TRASL 0
#define AVB 1
#define ADD 2
#define REMOVE 3
#define AVB_IN_IN 4

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
	vector *patches;
	matrix orient, orient_old;
	int index;

	int cell, cell_old;
	struct PatchyParticle *next;

	int n_aa;
	int type;
} PatchyParticle;

typedef struct {
	int N_side;
	int N;
	PatchyParticle **heads;
} LR_cells;

typedef struct {
	char configuration_prefix[512];
	char configuration_last[256];
	char sus_prefix[512];
	int print_csd;
	FILE *log;
	FILE *energy;
	FILE *density;
	FILE *acc;
} LR_IO;

typedef struct LR_system {
	int N, N_min, N_max;
	int substeps;
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
	int use_avb_in_in;
	int ensemble;
	void (*do_dynamics)(struct LR_system *, LR_IO *);

	llint *SUS_hist;
	llint **SUS_e_hist;
	double SUS_e_step;
	int SUS_e_bins;

	int overlap;

	int tries[N_MOVES];
	int accepted[N_MOVES];

	double disp_delta;
	double or_delta;

	double kf_delta_aa, kf_cosmax_aa, kf_sqr_rcut_aa;
	double avb_vin, avb_vout;
	double avb_p;

	LR_cells cells;

	int seed;
	PatchyParticle *particles;
} LR_system;

#endif /* DEFS_H_ */
