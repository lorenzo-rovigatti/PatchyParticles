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

void _init_cells(System *syst, Output *IO);

void rototraslate_particle(System *syst, PatchyParticle *p, vector disp, vector *orient);
void rollback_particle(System *syst, PatchyParticle *p);

void MC_add_remove(System *syst, Output *IO);
void MC_add_remove_single_cluster(System *syst, Output *IO);
void MC_avb_dyn(System *syst, Output *IO);
void MC_move_rototranslate(System *syst, Output *IO);

void init_MC(input_file *input, System *syst, Output *IO);
void make_initial_conf(System *syst, Output *IO, char *conf_name);
void fix_new_neighs(PatchyParticle *p);
double energy(System *syst, PatchyParticle *p);
void change_cell(System *syst, PatchyParticle *p);
void check_energy(System *syst, Output *IO);

#endif /* MC_H_ */
