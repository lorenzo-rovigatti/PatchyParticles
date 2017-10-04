/*
 * MC.c
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#include <float.h>
#include <math.h>

#include "MC.h"
#include "utils.h"
#include "LR_IO.h"
#include "neighs.h"

void init_MC(input_file *input, LR_system *syst, LR_IO *IO) {
	switch(syst->dynamics) {
	case 0:
		syst->do_dynamics = &MC_rototranslate_dyn;
		break;
	case 2:
		syst->do_dynamics = &MC_avb_dyn;
		break;
	default:
		die(IO, "Dynamics %d not supported\n", syst->dynamics);
		break;
	}
}

void rollback_particle(LR_system *syst, PatchyParticle *p) {
	p->r[0] = p->r_old[0];
	p->r[1] = p->r_old[1];
	p->r[2] = p->r_old[2];

	int i;
	for(i = 0; i < 3; i++) memcpy(p->orient[i], p->orient_old[i], sizeof(double)*3);
	for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(p->orient, syst->base_patches[i], p->patches[i]);

	// bring back the particle in the old cell
	if(p->cell != p->cell_old) {
		if(syst->cells.heads[p->cell]->index == p->index) syst->cells.heads[p->cell] = p->next;
		else {
			PatchyParticle *q = syst->cells.heads[p->cell];
			while(q->next != p) q = q->next;
			q->next = p->next;
		}

		PatchyParticle *old = syst->cells.heads[p->cell_old];
		syst->cells.heads[p->cell_old] = p;
		p->next = old;
		int c_old = p->cell;
		p->cell = p->cell_old;
		p->cell_old = c_old;
	}
}

void change_cell(LR_system *syst, PatchyParticle *p) {
	int ind[3];
	ind[0] = (int) ((p->r[0] / syst->L - floor(p->r[0] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[1] = (int) ((p->r[1] / syst->L - floor(p->r[1] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[2] = (int) ((p->r[2] / syst->L - floor(p->r[2] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);

	int cell_index = (ind[0] * syst->cells.N_side + ind[1]) * syst->cells.N_side + ind[2];
	if(cell_index == p->cell) {
		p->cell_old = p->cell;
		return;
	}

	// remove the particle from the old cell
	PatchyParticle *previous = NULL;
	PatchyParticle *current = syst->cells.heads[p->cell];
	assert(current != NULL);
	while(current->index != p->index) {
		previous = current;
		current = current->next;
		assert(previous->next->index == current->index);
	}
	if(previous == NULL) syst->cells.heads[p->cell] = p->next;
	else previous->next = p->next;

	// add the particle to the new cell
	p->next = syst->cells.heads[cell_index];
	syst->cells.heads[cell_index] = p;
	p->cell_old = p->cell;
	p->cell = cell_index;
}

void translate(LR_system *syst, PatchyParticle *p, vector disp) {
	p->r_old[0] = p->r[0];
	p->r_old[1] = p->r[1];
	p->r_old[2] = p->r[2];

	p->r[0] += disp[0];
	p->r[1] += disp[1];
	p->r[2] += disp[2];

	change_cell(syst, p);
}

void rototraslate(LR_system *syst, PatchyParticle *p, vector disp, vector *orient) {
	p->r_old[0] = p->r[0];
	p->r_old[1] = p->r[1];
	p->r_old[2] = p->r[2];

	p->r[0] += disp[0];
	p->r[1] += disp[1];
	p->r[2] += disp[2];

	int i, j;
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			p->orient_old[i][j] = p->orient[i][j];
			p->orient[i][j] = orient[i][j];
		}
	}

	for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(p->orient, syst->base_patches[i], p->patches[i]);

	change_cell(syst, p);
}

double current_energy(LR_system *syst, PatchyParticle *p) {
	syst->overlap = 0;
	return energy(syst, p);
}

double energy(LR_system *syst, PatchyParticle *p) {
	double E = 0.;
	if(syst->overlap) return E;

	int ind[3], loop_ind[3];
	ind[0] = (int) ((p->r[0] / syst->L - floor(p->r[0] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[1] = (int) ((p->r[1] / syst->L - floor(p->r[1] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[2] = (int) ((p->r[2] / syst->L - floor(p->r[2] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);

	int j, k, l, p_patch, q_patch;

	for(j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + syst->cells.N_side) % syst->cells.N_side;
		for(k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + syst->cells.N_side) % syst->cells.N_side;
			for(l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + syst->cells.N_side) % syst->cells.N_side;
				int loop_index = (loop_ind[0] * syst->cells.N_side + loop_ind[1]) * syst->cells.N_side + loop_ind[2];

				PatchyParticle *q = syst->cells.heads[loop_index];
				while(q != NULL) {
					if(q->index != p->index) {
						int val = kf_interact(syst, p, q, &p_patch, &q_patch);

						if(val == PATCH_BOND) {
							E -= 1.;
						}
						else if(val == -1) {
							vector dist = {q->r[0] - p->r[0], q->r[1] - p->r[1], q->r[2] - p->r[2]};
							dist[0] -= syst->L * rint(dist[0] / syst->L);
							dist[1] -= syst->L * rint(dist[1] / syst->L);
							dist[2] -= syst->L * rint(dist[2] / syst->L);

							double dist2 = SCALAR(dist, dist);
							if(dist2 > 0.999) syst->overlap = 2;
							else syst->overlap = 1;
							return 0.;
						}
					}
					q = q->next;
				}
			}
		}
	}

	return E;
}

void make_initial_conf(LR_system *syst, LR_IO *IO, char *conf_name) {
	int inserted = 0;
	while(inserted < syst->N) {
		// extract new position
		vector r = {drand48() * syst->L, drand48() * syst->L, drand48() * syst->L};
		PatchyParticle *p = syst->particles + inserted;
		p->r[0] = p->r[1] = p->r[2] = 0.;
		p->index = inserted;

		if(!would_overlap(syst, p, r)) {
			random_orientation(syst, p->orient);
			p->r[0] = r[0];
			p->r[1] = r[1];
			p->r[2] = r[2];

			// add the particle to the new cell
			int ind[3];
			ind[0] = (int) ((p->r[0] / syst->L - floor(p->r[0] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
			ind[1] = (int) ((p->r[1] / syst->L - floor(p->r[1] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
			ind[2] = (int) ((p->r[2] / syst->L - floor(p->r[2] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
			int cell_index = (ind[0] * syst->cells.N_side + ind[1]) * syst->cells.N_side + ind[2];
			p->next = syst->cells.heads[cell_index];
			syst->cells.heads[cell_index] = p;
			p->cell = p->cell_old = cell_index;

			inserted++;
		}
	}

	_print_conf(IO, syst, 0, conf_name);
}

void MC_avb_dyn(LR_system *syst, LR_IO *IO) {
	/*
	if(drand48() < syst->avb_p || syst->N < 2) MC_rototranslate_dyn(syst, IO);
	else {
		syst->tries[AVB]++;
		PatchyParticle *rec = syst->particles + (int)(drand48()*syst->N);
		int rec_n_neighs = rec->n_aa;
		int i;
		PatchyParticle *p;

		// bring in
		if(drand48() > 0.5) {
			// if all the particles but rec are in rec's bonding volume then no move is possible
			if(rec_n_neighs == (syst->N-1)) return;
			int found = 0;
			// select a particle which is not in rec's bonding volume
			do {
				p = syst->particles + (int)(drand48()*syst->N);
				found = 1;
				for(i = 0; i < syst->n_patches && found; i++) if(rec->neighs[i].p == p && rec->neighs[i].type == PATCH_BOND) found = 0;
			} while(!found || p == rec);

			vector new_r;
			matrix new_orient;
			// target patch
			int tarp = (int)(drand48() * syst->n_patches);
			place_inside_vbonding(syst, rec, new_r, new_orient, tarp);
			vector disp = {new_r[0] - p->r[0], new_r[1] - p->r[1], new_r[2] - p->r[2]};

			double deltaE = -current_energy(syst, p);
			rototraslate(syst, p, disp, new_orient);
			deltaE += energy(syst, p);

#ifdef DEBUG
			int p_patch = 0, q_patch = 0;
			assert(kf_interact(syst, rec, p, &p_patch, &q_patch) == PATCH_BOND);
#endif

			double acc = exp(-deltaE/syst->T) * (syst->N - rec_n_neighs - 1.) * syst->avb_vin / ((rec_n_neighs + 1.) * syst->avb_vout);
			if(!syst->overlap && drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;

				update_patch_bonds(syst, p);
			}
			else {
				rollback_particle(syst, p);
				syst->overlap = 0;
			}
		}
		// take out
		else {
			// we need rec to have neighbors in order to take one of them out...
			if(rec_n_neighs == 0) return;

			// pick a random particle from rec's neighbors
			int ip = (int) (drand48() * rec_n_neighs);
			int pc = -1;
			while(ip >= 0) {
				pc++;
				if(is_AA_bond(rec, pc)) ip--;
			}
			p = rec->neighs[pc].p;

#ifdef DEBUG
			int p_patch = 0, q_patch = 0;
			assert(p != NO_NEIGH);
			assert(kf_interact(syst, rec, p, &p_patch, &q_patch) == PATCH_BOND);
#endif

			vector new_r;
			matrix new_orient;
			vector new_patches[syst->n_patches];
			int buff;
			do {
				new_r[0] = drand48() * syst->L;
				new_r[1] = drand48() * syst->L;
				new_r[2] = drand48() * syst->L;
				random_orientation(syst, new_orient);
				for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(new_orient, syst->base_patches[i], new_patches[i]);
			} while(kf_would_interact(syst, rec, new_r, new_patches, &buff, &buff) == PATCH_BOND);

			vector disp = {new_r[0] - p->r[0], new_r[1] - p->r[1], new_r[2] - p->r[2]};
			if(would_overlap(syst, p, disp)) return;

			double deltaE = -current_energy(syst, p);
			rototraslate(syst, p, disp, new_orient);
			deltaE += energy(syst, p);

			double acc = exp(-deltaE/syst->T) * rec_n_neighs * syst->avb_vout / ((syst->N - rec_n_neighs) * syst->avb_vin);
			if(drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;

				update_patch_bonds(syst, p);
			}
			else {
				rollback_particle(syst, p);
				syst->overlap = 0;
			}
		}
	}
	*/
}

void MC_rototranslate_dyn(LR_system *syst, LR_IO *IO) {
	PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);
	int type = ROTO_TRASL;
	syst->tries[type]++;

	vector disp = {(drand48() - 0.5) * syst->disp_delta, (drand48() - 0.5) * syst->disp_delta, (drand48() - 0.5) * syst->disp_delta};
	// new orientation
	double t = drand48() * syst->or_delta;
	vector axis;
	matrix new_orient;
	random_vector_on_sphere(axis);
	rotate_orient(p->orient, new_orient, axis, t);

	// apply changes to p
	double deltaE = -current_energy(syst, p);
	rototraslate(syst, p, disp, new_orient);
	deltaE += energy(syst, p);

	if(!syst->overlap && (deltaE < 0. || drand48() < exp(-deltaE/syst->T))) {
		syst->energy += deltaE;
		syst->accepted[type]++;
	}
	else {
		rollback_particle(syst, p);
		syst->overlap = 0;
	}
}

void MC_add_remove(LR_system *syst, LR_IO *IO) {
	// try to add a particle
	if(drand48() < 0.5) {
		if(syst->N == syst->N_max) return;
		syst->tries[ADD]++;

		PatchyParticle *p = syst->particles + syst->N;
		p->index = syst->N;

		double acc_factor = 1;
		if(syst->N == 0) {
			p->r[0] = drand48() * syst->L;
			p->r[1] = drand48() * syst->L;
			p->r[2] = drand48() * syst->L;

			random_orientation(syst, p->orient);
			acc_factor = syst->V / (syst->N + 1.);
		}
		// or put the new particle in the bonding volume of an
		// already present particle
		else {
			PatchyParticle *rec = syst->particles + (int)(drand48()*syst->N);
			int tarp = (int)(drand48() * syst->n_patches);
			place_inside_vbonding(syst, rec, p->r, p->orient, tarp);
			acc_factor = syst->N * syst->avb_vin / (syst->N + 1.);
#ifdef DEBUG
			int p_patch = 0, q_patch = 0, i;
			for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(p->orient, syst->base_patches[i], p->patches[i]);
			int var = kf_interact(syst, rec, p, &p_patch, &q_patch);
			if(var != PATCH_BOND) printf("%d\n", var);
			assert(kf_interact(syst, rec, p, &p_patch, &q_patch) == PATCH_BOND);
#endif
		}

		int i;
		for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(p->orient, syst->base_patches[i], p->patches[i]);

		double delta_E = energy(syst, p);
		double acc = exp(-delta_E/syst->T) * syst->z * acc_factor;

		if(!syst->overlap && drand48() < acc) {
			syst->energy += delta_E;
			p->n_aa = 0;
			p->type = TYPE_NO;

			// add the particle to the new cell
			int ind[3];
			ind[0] = (int) ((p->r[0] / syst->L - floor(p->r[0] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
			ind[1] = (int) ((p->r[1] / syst->L - floor(p->r[1] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
			ind[2] = (int) ((p->r[2] / syst->L - floor(p->r[2] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
			int cell_index = (ind[0] * syst->cells.N_side + ind[1]) * syst->cells.N_side + ind[2];
			p->next = syst->cells.heads[cell_index];
			syst->cells.heads[cell_index] = p;
			p->cell = p->cell_old = cell_index;

			syst->N++;
			syst->accepted[ADD]++;
		}
	}
	// try to remove a particle
	else {
		if(syst->N == syst->N_min) return;
		syst->tries[REMOVE]++;

		PatchyParticle *p = syst->particles + (int)(drand48() * syst->N);

		double delta_E = -current_energy(syst, p);

		double acc_factor = syst->N / syst->V;
		double acc = exp(-delta_E/syst->T) * acc_factor / syst->z;
		if(drand48() < acc) {
			syst->energy += delta_E;
			syst->N--;

			// remove the particle from the old cell
			PatchyParticle *previous = NULL;
			PatchyParticle *current = syst->cells.heads[p->cell];
			assert(current != NULL);
			while(current->index != p->index) {
				previous = current;
				current = current->next;
				assert(previous->next->index == current->index);
			}
			if(previous == NULL) syst->cells.heads[p->cell] = p->next;
			else previous->next = p->next;

			if(p->index != syst->N) {
				PatchyParticle *q = syst->particles + syst->N;
				assert(q->index != p->index);
				// we have to change the last particle's identity
				// first we remove it from its cell
				previous = NULL;
				current = syst->cells.heads[q->cell];
				assert(current != NULL);
				while(current->index != q->index) {
					previous = current;
					current = current->next;
					assert(previous->next->index == current->index);
				}
				if(previous == NULL) syst->cells.heads[q->cell] = q->next;
				else previous->next = q->next;

				// copy its type, position and patches onto p's memory position
				p->r[0] = q->r[0];
				p->r[1] = q->r[1];
				p->r[2] = q->r[2];

				p->n_aa = q->n_aa;
				p->type = q->type;

				int i;
				for(i = 0; i < 3; i++) memcpy(p->orient[i], q->orient[i], sizeof(double) * 3);
				for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(p->orient, syst->base_patches[i], p->patches[i]);

				// and finally add it back to its former cell
				p->cell = q->cell;
				p->next = syst->cells.heads[p->cell];
				syst->cells.heads[p->cell] = p;
			}

			syst->accepted[REMOVE]++;
		}
	}
}

void MC_step(LR_system *syst, LR_IO *IO) {
	int i;
	for(i = 0; i < syst->substeps; i++) {
		syst->overlap = 0;
		// TODO: use function pointers to incapsulate the ensembles
		if(syst->ensemble != 0 && drand48() < 0.01) {
			MC_add_remove(syst, IO);
			syst->substeps = syst->N*10;
			if(syst->substeps < 1000) syst->substeps = 1000;
			if(syst->ensemble == 3) syst->SUS_hist[syst->N - syst->N_min]++;
			// TODO: remove the ensemble == 6 (aka the 2D histogram)
			else if(syst->ensemble == 6) {
				int epos = (int)(fabs(syst->energy) / syst->SUS_e_step);
				syst->SUS_e_hist[syst->N - syst->N_min][epos]++;
			}
		}
		else if(syst->N > 0) syst->do_dynamics(syst, IO);
	}
}

void check_energy(LR_system *syst, LR_IO *IO) {
	int i;
	double E = 0.;
	syst->overlap = 0;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		E += current_energy(syst, p);
		gram_schmidt(p->orient[0], p->orient[1], p->orient[2]);
	}
	E *= 0.5;

	if(syst->overlap) die(IO, "\nComputing energy from scratch resulted in an overlap, aborting");
	if(fabs(syst->energy) > 1e-5 && fabs((E - syst->energy)/syst->energy) > 1e-5) {
		die(IO, "\nEnergy check failed, old energy = %lf, new energy = %lf\n", syst->energy, E);
	}

	syst->energy = E;
}
