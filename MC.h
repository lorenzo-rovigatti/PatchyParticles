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

void init_cells(LR_system *syst, LR_IO *IO);

void rototraslate_particle(LR_system *syst, PatchyParticle *p, vector disp, vector *orient);
void rollback_particle(LR_system *syst, PatchyParticle *p);

void MC_add_remove(LR_system *syst, LR_IO *IO);
void MC_add_remove_single_cluster(LR_system *syst, LR_IO *IO);
void MC_avb_dyn(LR_system *syst, LR_IO *IO);
void MC_move_rototranslate(LR_system *syst, LR_IO *IO);

void init_MC(input_file *input, LR_system *syst, LR_IO *IO);
void make_initial_conf(LR_system *syst, LR_IO *IO, char *conf_name);
void MC_sweep(LR_system *syst, LR_IO *IO);
void fix_new_neighs(PatchyParticle *p);
double energy(LR_system *syst, PatchyParticle *p);
void change_cell(LR_system *syst, PatchyParticle *p);
void check_energy(LR_system *syst, LR_IO *IO);

#endif /* MC_H_ */
