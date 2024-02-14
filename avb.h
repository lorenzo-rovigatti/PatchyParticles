#ifndef AVB_H_
#define AVB_H_

typedef struct Output Output;
typedef struct System System;
typedef struct PatchyParticle PatchyParticle;
typedef struct input_file input_file;

typedef struct _avbmc {
	PatchyParticle **neighbours;
	int *num_neighbours_patch;
	int num_neighbours;

	double avb_vin, avb_vout;
	double avb_vin_per_unit;
	double avb_p;
} avbmc;

void AVBMC_init(input_file *input, System *syst, Output *IO);
void AVBMC_free();

void AVBMC_dynamics(System *syst, Output *IO);
void AVBMC_dynamics_colors(System *syst, Output *IO);

#endif
