#ifndef VMMC_H_
#define VMMC_H_

#define VMMC_ROTATION (1)
#define VMMC_TRANSLATION (2)

#include "defs.h"

typedef struct input_file input_file;

typedef struct vmmc_d {
	double max_move;
	int max_cluster;
	PatchyParticle ** possible_links;
	int n_possible_links;
	PatchyParticle ** prelinked_particles;
	int n_prelinked_particles;
	PatchyParticle ** clust;
	int n_clust;
	int * is_in_cluster;
	int which_move;
	matrix rotation;
} vmmc_d;

void VMMC_init(input_file *input, System *syst, Output *IO);
void VMMC_free();

void VMMC_dynamics(System *syst, Output *IO);

#endif
