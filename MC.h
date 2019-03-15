/*
 * MC.h
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#ifndef MC_H_
#define MC_H_

#include "defs.h"

typedef struct input_file input_file;

void MC_rototraslate_particle(System *syst, PatchyParticle *p, vector disp, vector *orient);
void MC_rollback_particle(System *syst, PatchyParticle *p);

int MC_would_interact(System *syst, PatchyParticle *p, vector r, matrix orient, int *onp, int *onq);
int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, int *onp, int *onq);

double MC_energy(System *syst, PatchyParticle *p);

void MC_add_remove(System *syst, Output *output_files);
void MC_move_rototranslate(System *syst, Output *output_files);
void MC_change_volume(System *syst, Output *output_files);
void MC_change_Lx(System *syst, Output *output_files);

void MC_init(input_file *input, System *syst, Output *output_files);
void MC_add_remove_biased(System *syst, Output *IO);
void MC_add_remove_single_cluster(System *syst, Output *IO);
void MC_move_rototranslate(System *syst, Output *IO);

void bsus_update_histo(System *syst);

void MC_init(input_file *input, System *syst, Output *IO);
void MC_free(System *syst);
void MC_change_cell(System *syst, PatchyParticle *p);
void MC_check_energy(System *syst, Output *output_files);

#endif /* MC_H_ */
