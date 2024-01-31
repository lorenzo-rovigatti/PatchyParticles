/*
 * defs.h
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#ifndef DEFS_H_
#define DEFS_H_

#define N_MOVES 10
#define ROTO_TRASL 0
#define AVB 1
#define ADD 2
#define REMOVE 3
#define AVB_IN_IN 4
#define MOVE_VMMC 5
#define VOLUME 6
#define LX 7
#define TRANSFER 8
#define USCNTMOVE 9

#define NVT 0
#define GC 1
#define SUS 3
#define NPT 5
#define BSUS 4
#define GIBBS 6
#define CNTUS 7

#define RTMC 0
#define VMMC 1
#define AVBMC 2

#define OVERLAP -1
#define NO_BOND 0
#define PATCH_BOND 1

#include "cells.h"

#include <stdio.h>
#include <complex.h>
#include <assert.h>
#include <signal.h>

typedef double vector[3];
typedef double matrix[3][3];
typedef long long int llint;
typedef struct Output Output;

typedef struct PatchyParticle {
	vector r, r_old;
	matrix orientation, orientation_old;
	int index;

	int n_patches;
	vector *patches, *base_patches;

	int cell, cell_old;
	int specie;
} PatchyParticle;

typedef struct System {
	int N, N_min, N_max;
	vector box;
	double V;
	double T;
	double z;
	double P;
	double energy;
	int n_patches;

	/**
	 * This matrix is initialised so as to transform base_patches[0] to the 0, 0, 1 vector
	 */
	matrix base_orient;
	vector *base_patches;

	int dynamics;
	int ensemble;
	int Lx_move;
	double yz_ratio;
	int piston_direction;
	void (*do_dynamics)(struct System *, Output *);
	void (*do_ensemble)(struct System *, Output *);

	// non-biased sus variables
	llint *SUS_hist;
	llint **SUS_e_hist;
	double SUS_e_step;
	int SUS_e_bins;

	// biased sus variables
	double *bsus_collect;
	double *bsus_tm;
	double *bsus_normvec;
	double *bsus_pm;


	int overlap;

	int tries[N_MOVES];
	int accepted[N_MOVES];

	double disp_max;
	double theta_max;
	double rescale_factor_max;
	double Lx_change_max;
	double Lyz_min;
	double Lyz_max;


	double kf_delta, kf_cosmax, kf_sqr_rcut;
	double r_cut;

	Cells *cells;
	Cells *old_cells;

	int seed;
	PatchyParticle *particles;

	// COLORS ///////////////////////
	int num_species;
	int num_colors;
	int **colorint;
	int *ncolorint;
	int **particlescolor; // particlescolor[specie,patch]=color
	int **color; // color[color,specie]=how many patches of that color
	int *species_count; // species_count[specie]=how many particles of that specie
	int **bonding_volume_units;
	/////////////////////////////////

	// GIBBS  ///////////////////////
	double gibbsVolumeFrequency;
	double gibbsSwapFrequency;
	double gibbsVolumeDeltamax;
	/////////////////////////////////

	// UMBRELLA SAMPLING
	double US_OP;
	double US_OP_MIN;
	double US_OP_MAX;
	double *US_old_pos;
	double *US_old_orientation;
	int US_bias_type;
	double US_k;
	int US_steps;
	double US_OP_0;
		

} System;

#endif /* DEFS_H_ */
