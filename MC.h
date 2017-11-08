/*
 * MC.h
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#ifndef MC_H_
#define MC_H_

#include "defs.h"
#include "parse_input.h"

void rototraslate_particle(System *syst, PatchyParticle *p, vector disp, vector *orient);
void rollback_particle(System *syst, PatchyParticle *p);

int MC_would_interact(System *syst, PatchyParticle *p, vector r, matrix orient, int *onp, int *onq);
int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, int *onp, int *onq);

void MC_add_remove(System *syst, Output *IO);
void MC_change_volume(System *syst, Output *IO);
void MC_move_rototranslate(System *syst, Output *IO);

void MC_init(input_file *input, System *syst, Output *IO);
void MC_free(System *syst);
double energy(System *syst, PatchyParticle *p);
void change_cell(System *syst, PatchyParticle *p);
void check_energy(System *syst, Output *IO);

#endif /* MC_H_ */
