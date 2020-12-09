/*
 * MC.c
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#include "avb.h"
#include "MC.h"
#include "output.h"
#include "utils.h"
#include "vmmc.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void do_NVT(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N; i++) {
		syst->do_dynamics(syst, output_files);

		double R = drand48();
		if(syst->Lx_move && R < 1. / syst->N) {
			MC_change_Lx(syst, output_files);
		}
	}
}

void do_GC(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N_max; i++) {
		double R = drand48();
		if(R < 0.01) {
			MC_add_remove(syst, output_files);
		}
		else if(syst->Lx_move && R < (0.01 + 1. / syst->N)) {
			MC_change_Lx(syst, output_files);
		}
		else if(syst->N > 0) syst->do_dynamics(syst, output_files);
	}
}

void do_SUS(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N_max; i++) {
		double R = drand48();
		if(R < 0.01) {
			MC_add_remove(syst, output_files);
			syst->SUS_hist[syst->N - syst->N_min]++;
		}
		else if(syst->Lx_move && R < (0.01 + 1. / syst->N)) {
			MC_change_Lx(syst, output_files);
		}
		else if(syst->N > 0) syst->do_dynamics(syst, output_files);
	}
}

void do_NPT(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N; i++) {
		if(drand48() < 1. / syst->N) {
			MC_change_volume(syst, output_files);
		}
		else if(syst->N > 0) syst->do_dynamics(syst, output_files);
	}
}


void do_BSUS(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N_max; i++) {
		double R = drand48();
		if(R < 0.5) {
			MC_add_remove_biased(syst, output_files);
		}
		else if(syst->Lx_move && R < (0.5 + 1. / syst->N)) {
			MC_change_Lx(syst, output_files);
		}
		else if(syst->N > 0) syst->do_dynamics(syst, output_files);
	}
}

void MC_init(input_file *input, System *syst, Output *IO) {
	/**
	 * Here we set the pointer to the function that will be used to make a Monte Carlo step
	 * according to the ensemble specified in the input file
	 */
	switch(syst->ensemble) {
	case NVT:
		syst->do_ensemble = &do_NVT;
		break;
	case GC:
		syst->do_ensemble = &do_GC;
		break;
	case SUS:
		syst->do_ensemble = &do_SUS;
		break;
	case NPT:
		syst->do_ensemble = &do_NPT;
		break;
	case BSUS:
		syst->do_ensemble= &do_BSUS;
		break;
	default:
		output_exit(IO, "Ensemble %d not supported\n", syst->ensemble);
		break;
	}

	/**
	 * Here we set the pointer to the function that will be used to perform the type of dynamics
	 * specified in the input file
	 */
	switch(syst->dynamics) {
	case RTMC:
		syst->do_dynamics = &MC_move_rototranslate;
		break;
	case VMMC:
		VMMC_init(input, syst, IO);
		syst->do_dynamics = &VMMC_dynamics;
		break;
	case AVBMC:
		AVBMC_init(input, syst, IO);
		syst->do_dynamics = &AVBMC_dynamics;
		break;
	default:
		output_exit(IO, "Dynamics %d not supported\n", syst->dynamics);
		break;
	}
}

void MC_free(System *syst) {
	switch(syst->dynamics) {
	case RTMC:
		break;
	case VMMC:
		VMMC_free();
		break;
	case AVBMC:
		AVBMC_free();
		break;
	}
}

void MC_rollback_particle(System *syst, PatchyParticle *p) {
	p->r[0] = p->r_old[0];
	p->r[1] = p->r_old[1];
	p->r[2] = p->r_old[2];

	int i;
	for(i = 0; i < 3; i++) {
		memcpy(p->orientation[i], p->orientation_old[i], sizeof(double) * 3);
	}
	for(i = 0; i < p->n_patches; i++) {
		MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
	}

	Cells *cells = syst->cells;
	// bring the particle back in the old cell
	if(p->cell != p->cell_old) {
		if(cells->heads[p->cell]->index == p->index) cells->heads[p->cell] = cells->next[p->index];
		else {
			PatchyParticle *q = cells->heads[p->cell];
			while(cells->next[q->index] != p) {
				q = cells->next[q->index];
			}
			cells->next[q->index] = cells->next[p->index];
		}

		PatchyParticle *old = cells->heads[p->cell_old];
		cells->heads[p->cell_old] = p;
		cells->next[p->index] = old;
		int c_old = p->cell;
		p->cell = p->cell_old;
		p->cell_old = c_old;
	}
}

void MC_change_cell(System *syst, PatchyParticle *p) {
	int ind[3];
	int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
	if(cell_index == p->cell) {
		p->cell_old = p->cell;
		return;
	}

	// remove the particle from the old cell
	Cells *cells = syst->cells;
	PatchyParticle *previous = NULL;
	PatchyParticle *current = cells->heads[p->cell];
	assert(current != NULL);
	while(current->index != p->index) {
		previous = current;
		current = cells->next[current->index];
		assert(cells->next[previous->index]->index == current->index);
	}
	if(previous == NULL) cells->heads[p->cell] = cells->next[p->index];
	else cells->next[previous->index] = cells->next[p->index];

	// add the particle to the new cell
	cells->next[p->index] = cells->heads[cell_index];
	cells->heads[cell_index] = p;
	p->cell_old = p->cell;
	p->cell = cell_index;
}

void MC_rototraslate_particle(System *syst, PatchyParticle *p, vector disp, vector *orient) {
	p->r_old[0] = p->r[0];
	p->r_old[1] = p->r[1];
	p->r_old[2] = p->r[2];

	p->r[0] += disp[0];
	p->r[1] += disp[1];
	p->r[2] += disp[2];

	int i, j;
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			p->orientation_old[i][j] = p->orientation[i][j];
			p->orientation[i][j] = orient[i][j];
		}
	}

	for(i = 0; i < p->n_patches; i++) {
		MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
	}

	MC_change_cell(syst, p);
}

int MC_would_interact(System *syst, PatchyParticle *p, vector r, vector *patches, int *onp, int *onq) {
	vector dist = {r[0] - p->r[0], r[1] - p->r[1], r[2] - p->r[2]};
	dist[0] -= syst->box[0] * rint(dist[0] / syst->box[0]);
	dist[1] -= syst->box[1] * rint(dist[1] / syst->box[1]);
	dist[2] -= syst->box[2] * rint(dist[2] / syst->box[2]);

	double dist2 = SCALAR(dist, dist);

	if(dist2 < 1.) return OVERLAP;
	else if(dist2 < syst->kf_sqr_rcut) {
		// versor
		double norm = sqrt(dist2);
		dist[0] /= norm;
		dist[1] /= norm;
		dist[2] /= norm;

		int pp, pq;
		for(pp = 0; pp < syst->n_patches; pp++) {
			double p_cos = SCALAR(dist, p->patches[pp]);

			if(p_cos > syst->kf_cosmax) {
				for(pq = 0; pq < syst->n_patches; pq++) {
					double q_cos = -SCALAR(dist, patches[pq]);
					if(q_cos > syst->kf_cosmax) {
						*onp = pp;
						*onq = pq;
						return PATCH_BOND;
					}
				}
			}
		}
	}

	return NO_BOND;
}

inline int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, int *onp, int *onq) {
	return MC_would_interact(syst, p, q->r, q->patches, onp, onq);
}

double MC_energy(System *syst, PatchyParticle *p) {
	syst->overlap = 0;
	double E = 0.;

	int ind[3], loop_ind[3];
	cells_fill_and_get_idx_from_particle(syst, p, ind);

	int j, k, l, p_patch, q_patch;

	for(j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + syst->cells->N_side[0]) % syst->cells->N_side[0];
		for(k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + syst->cells->N_side[1]) % syst->cells->N_side[1];
			for(l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + syst->cells->N_side[2]) % syst->cells->N_side[2];
				int loop_index = (loop_ind[0] * syst->cells->N_side[1] + loop_ind[1]) * syst->cells->N_side[2] + loop_ind[2];

				PatchyParticle *q = syst->cells->heads[loop_index];
				while(q != NULL) {
					if(q->index != p->index) {
						int val = MC_interact(syst, p, q, &p_patch, &q_patch);

						if(val == PATCH_BOND) {
							E -= 1.;
						}
						else if(val == OVERLAP) {
							syst->overlap = 1;
							return 0.;
						}
					}
					q = syst->cells->next[q->index];
				}
			}
		}
	}

	return E;
}

void MC_move_rototranslate(System *syst, Output *IO) {
	PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);
	int type = ROTO_TRASL;
	syst->tries[type]++;

	vector disp;
	disp[0] = (drand48() - 0.5) * syst->disp_max;
	disp[1] = (drand48() - 0.5) * syst->disp_max;
	disp[2] = (drand48() - 0.5) * syst->disp_max;
	// new orientation
	double theta = drand48() * syst->theta_max;
	vector axis;
	matrix new_orient;
	random_vector_on_sphere(axis);
	utils_rotate_matrix(p->orientation, new_orient, axis, theta);

	// apply changes to p
	double deltaE = -MC_energy(syst, p);
	MC_rototraslate_particle(syst, p, disp, new_orient);
	deltaE += MC_energy(syst, p);

	if(!syst->overlap && (deltaE < 0. || drand48() < exp(-deltaE / syst->T))) {
		syst->energy += deltaE;
		syst->accepted[type]++;
	}
	else {
		MC_rollback_particle(syst, p);
		syst->overlap = 0;
	}
}

void MC_add_remove(System *syst, Output *IO) {
	// try to add a particle
	if(drand48() < 0.5) {
		if(syst->N == syst->N_max) {
			if(syst->ensemble == GC) output_exit(IO, "The system contains the maximum number of particles set in the input file (%d), increase GC_N_max\n", syst->N_max);
			return;
		}
		syst->tries[ADD]++;

		PatchyParticle *p = syst->particles + syst->N;
		p->index = syst->N;

		p->r[0] = drand48() * syst->box[0];
		p->r[1] = drand48() * syst->box[1];
		p->r[2] = drand48() * syst->box[2];

		random_orientation(syst, p->orientation);

		int i;
		for(i = 0; i < p->n_patches; i++) {
			MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
		}

		double delta_E = MC_energy(syst, p);
		double acc = exp(-delta_E / syst->T) * syst->z * syst->V / (syst->N + 1.);

		if(!syst->overlap && drand48() < acc) {
			syst->energy += delta_E;

			// add the particle to the new cell
			int ind[3];
			int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
			syst->cells->next[p->index] = syst->cells->heads[cell_index];
			syst->cells->heads[cell_index] = p;
			p->cell = p->cell_old = cell_index;

			syst->N++;
			syst->accepted[ADD]++;
		}
		else {
			syst->overlap = 0;
		}
	}
	// try to remove a particle
	else {
		if(syst->N == syst->N_min) return;
		syst->tries[REMOVE]++;

		PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);

		double delta_E = -MC_energy(syst, p);
		double acc = exp(-delta_E / syst->T) * syst->N / (syst->V * syst->z);
		if(drand48() < acc) {
			syst->energy += delta_E;
			syst->N--;

			// remove the particle from the old cell
			PatchyParticle *previous = NULL;
			PatchyParticle *current = syst->cells->heads[p->cell];
			assert(current != NULL);
			while(current->index != p->index) {
				previous = current;
				current = syst->cells->next[current->index];
				assert(syst->cells->next[previous->index]->index == current->index);
			}
			if(previous == NULL) syst->cells->heads[p->cell] = syst->cells->next[p->index];
			else syst->cells->next[previous->index] = syst->cells->next[p->index];

			if(p->index != syst->N) {
				PatchyParticle *q = syst->particles + syst->N;
				assert(q->index != p->index);
				// we have to change the last particle's identity
				// first we remove it from its cell
				previous = NULL;
				current = syst->cells->heads[q->cell];
				assert(current != NULL);
				while(current->index != q->index) {
					previous = current;
					current = syst->cells->next[current->index];
					assert(syst->cells->next[previous->index]->index == current->index);
				}
				if(previous == NULL) syst->cells->heads[q->cell] = syst->cells->next[q->index];
				else syst->cells->next[previous->index] = syst->cells->next[q->index];

				// copy its type, position and patches onto p's memory position
				p->r[0] = q->r[0];
				p->r[1] = q->r[1];
				p->r[2] = q->r[2];

				int i;
				for(i = 0; i < 3; i++) {
					memcpy(p->orientation[i], q->orientation[i], sizeof(double) * 3);
				}
				for(i = 0; i < p->n_patches; i++) {
					MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
				}

				// and finally add it back to its former cell
				p->cell = q->cell;
				syst->cells->next[p->index] = syst->cells->heads[p->cell];
				syst->cells->heads[p->cell] = p;
			}

			syst->accepted[REMOVE]++;
		}
	}
}

void MC_change_volume(System *syst, Output *IO) {
	syst->tries[VOLUME]++;

	double delta_E = -syst->energy;
	double initial_V = syst->box[0] * syst->box[1] * syst->box[2];

	// pick a random direction
	int dir = drand48() * 3;
	// extract a random volume change
	double ln_final_V = log(initial_V) + (drand48() - 0.5)*syst->rescale_factor_max;
	double final_V = exp(ln_final_V);

	double old_side = syst->box[dir];
	double new_side = final_V;
	int i = 0;
	for(i = 0; i < 3; i++) {
		if(i != dir) {
			new_side /= syst->box[i];
		}
	}

	syst->box[dir] = new_side;
	double rescale_factor = new_side / old_side;

	// if we compress the system too much we'll have to recompute the cells
	if((syst->box[dir] / syst->cells->N_side[dir]) < syst->r_cut) {
		cells_free(syst->cells);
		cells_init(syst, IO, syst->r_cut);
		cells_fill(syst);
	}

	// rescale particles' positions
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		p->r[dir] *= rescale_factor;
	}

	// compute the new energy
	int overlap_found = 0;
	for(i = 0; i < syst->N && !overlap_found; i++) {
		PatchyParticle *p = syst->particles + i;
		delta_E += MC_energy(syst, p) * 0.5;
		if(syst->overlap) {
			overlap_found = 1;
		}
	}

	double delta_V = final_V - initial_V;
	double exp_arg = -(syst->P * delta_V + delta_E) / syst->T + (syst->N + 1) * log(rescale_factor);
	if(!overlap_found && drand48() < exp(exp_arg)) {
		syst->V = final_V;
		syst->energy += delta_E;
		syst->accepted[VOLUME]++;
	}
	else {
//		if(!overlap_found) printf("rejected %lf -- %lf %lf %lf\n", syst->P*delta_V, delta_E, (syst->N + 1) * log(rescale_factor), exp(exp_arg));
		for(i = 0; i < syst->N; i++) {
			PatchyParticle *p = syst->particles + i;
			p->r[dir] /= rescale_factor;
		}

		syst->box[dir] = old_side;
		syst->overlap = 0;
	}
}

void MC_change_Lx(System *syst, Output *IO) {
	syst->tries[LX]++;

	double delta_E = -syst->energy;

	vector old_sides;
	set_vector(old_sides, syst->box[0], syst->box[1], syst->box[2]);
	syst->box[0] += (drand48() - 0.5) * syst->Lx_change_max;
	// we keep the volume constant
	syst->box[1] = syst->box[2] = sqrt(syst->V / syst->box[0]);

	// early rejection if Ly (or, equivalently, Lz) is not within the allowed range
	if(syst->box[1] < syst->Lyz_min || syst->box[1] > syst->Lyz_max) {
		set_vector(syst->box, old_sides[0], old_sides[1], old_sides[2]);
		return;
	}

	// if we compress the system too much we'll have to recompute the cells
	if((syst->box[0] / syst->cells->N_side[0]) < syst->r_cut || (syst->box[1] / syst->cells->N_side[1]) < syst->r_cut) {
		cells_free(syst->cells);
		cells_init(syst, IO, syst->r_cut);
		cells_fill(syst);
	}

	vector rescale_factors = {
		syst->box[0] / old_sides[0],
		syst->box[1] / old_sides[1],
		syst->box[2] / old_sides[2]
	};

	// rescale particles' positions
	int i;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		p->r[0] *= rescale_factors[0];
		p->r[1] *= rescale_factors[1];
		p->r[2] *= rescale_factors[2];
	}

	// compute the new energy
	int overlap_found = 0;
	for(i = 0; i < syst->N && !overlap_found; i++) {
		PatchyParticle *p = syst->particles + i;
		delta_E += MC_energy(syst, p) * 0.5;
		if(syst->overlap) overlap_found = 1;
	}

	if(!overlap_found && drand48() < exp(delta_E / syst->T)) {
		syst->energy += delta_E;
		syst->accepted[LX]++;
	}
	else {
		for(i = 0; i < syst->N; i++) {
			PatchyParticle *p = syst->particles + i;
			p->r[0] /= rescale_factors[0];
			p->r[1] /= rescale_factors[1];
			p->r[2] /= rescale_factors[2];
		}

		set_vector(syst->box, old_sides[0], old_sides[1], old_sides[2]);

		// we might have to recompute the cells once again
		if((syst->box[0] / syst->cells->N_side[0]) < syst->r_cut || (syst->box[1] / syst->cells->N_side[1]) < syst->r_cut) {
			cells_free(syst->cells);
			cells_init(syst, IO, syst->r_cut);
			cells_fill(syst);
		}

		syst->overlap = 0;
	}
}

void MC_add_remove_biased(System *syst, Output *IO) {
	// try to add a particle

	if(drand48() < 0.5) {
		if(syst->N == syst->N_max) {
			if(syst->ensemble == GC) output_exit(IO, "The system contains the maximum number of particles set in the input file (%d), increase GC_N_max\n", syst->N_max);

			syst->bsus_collect[3*(syst->N-syst->N_min)+1]+=1.;

			return;
		}
		syst->tries[ADD]++;

		PatchyParticle *p = syst->particles + syst->N;
		p->index = syst->N;

		p->r[0] = drand48() * syst->box[0];
		p->r[1] = drand48() * syst->box[1];
		p->r[2] = drand48() * syst->box[2];

		random_orientation(syst, p->orientation);

		int i;
		for(i = 0; i < p->n_patches; i++) {
			MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
		}

		double delta_E = MC_energy(syst, p);
		double acc = exp(-delta_E / syst->T) * syst->z * syst->V / (syst->N + 1.);

		double pa;
		if (syst->overlap==1)
			pa=0;
		else
			pa=(acc>1 ? 1 : acc);

		syst->bsus_collect[3*(syst->N-syst->N_min)+2]+=pa;
		syst->bsus_collect[3*(syst->N-syst->N_min)+1]+=(1.-pa);

		double bias=syst->bsus_pm[syst->N-syst->N_min+1]-syst->bsus_pm[syst->N-syst->N_min];


		if(!syst->overlap && drand48() < exp(-bias)*acc) {
			syst->energy += delta_E;

			// add the particle to the new cell
			int ind[3];
			int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
			syst->cells->next[p->index] = syst->cells->heads[cell_index];
			syst->cells->heads[cell_index] = p;
			p->cell = p->cell_old = cell_index;

			syst->N++;
			syst->accepted[ADD]++;
		}
		else {
			syst->overlap = 0;
		}
	}
	// try to remove a particle
	else {
		if(syst->N == syst->N_min)
		{
			syst->bsus_collect[1]+=1.;
			return;
		}
		syst->tries[REMOVE]++;

		PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);

		double delta_E = -MC_energy(syst, p);
		double acc = exp(-delta_E / syst->T) * syst->N / (syst->V * syst->z);

		double pa=(acc>1 ? 1 : acc);

		syst->bsus_collect[3*(syst->N-syst->N_min)+0]+=pa;
		syst->bsus_collect[3*(syst->N-syst->N_min)+1]+=(1.-pa);

		double bias=syst->bsus_pm[syst->N-syst->N_min-1]-syst->bsus_pm[syst->N-syst->N_min];

		if(drand48() < exp(-bias)*acc) {
			syst->energy += delta_E;
			syst->N--;

			// remove the particle from the old cell
			PatchyParticle *previous = NULL;
			PatchyParticle *current = syst->cells->heads[p->cell];
			assert(current != NULL);
			while(current->index != p->index) {
				previous = current;
				current = syst->cells->next[current->index];
				assert(syst->cells->next[previous->index]->index == current->index);
			}
			if(previous == NULL) syst->cells->heads[p->cell] = syst->cells->next[p->index];
			else syst->cells->next[previous->index] = syst->cells->next[p->index];

			if(p->index != syst->N) {
				PatchyParticle *q = syst->particles + syst->N;
				assert(q->index != p->index);
				// we have to change the last particle's identity
				// first we remove it from its cell
				previous = NULL;
				current = syst->cells->heads[q->cell];
				assert(current != NULL);
				while(current->index != q->index) {
					previous = current;
					current = syst->cells->next[current->index];
					assert(syst->cells->next[previous->index]->index == current->index);
				}
				if(previous == NULL) syst->cells->heads[q->cell] = syst->cells->next[q->index];
				else syst->cells->next[previous->index] = syst->cells->next[q->index];

				// copy its type, position and patches onto p's memory position
				p->r[0] = q->r[0];
				p->r[1] = q->r[1];
				p->r[2] = q->r[2];

				int i;
				for(i = 0; i < 3; i++) {
					memcpy(p->orientation[i], q->orientation[i], sizeof(double) * 3);
				}
				for(i = 0; i < p->n_patches; i++) {
					MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
				}

				// and finally add it back to its former cell
				p->cell = q->cell;
				syst->cells->next[p->index] = syst->cells->heads[p->cell];
				syst->cells->heads[p->cell] = p;
			}

			syst->accepted[REMOVE]++;
		}
	}
}


void bsus_update_histo(System *syst)
{

	int histogram_size=syst->N_max-syst->N_min+1;


	syst->bsus_pm[0]=1.0;

	int i;

	for (i=0;i<histogram_size;i++)
	{
		syst->bsus_normvec[i]=syst->bsus_collect[3*i];
		syst->bsus_normvec[i]+=syst->bsus_collect[3*i+1];
		syst->bsus_normvec[i]+=syst->bsus_collect[3*i+2];
	}

	for (i=0;i<histogram_size;i++)
	{
		if (syst->bsus_normvec[i]>0)
		{
			syst->bsus_tm[3*i+0]=syst->bsus_collect[3*i+0]/syst->bsus_normvec[i];
			syst->bsus_tm[3*i+1]=syst->bsus_collect[3*i+1]/syst->bsus_normvec[i];
			syst->bsus_tm[3*i+2]=syst->bsus_collect[3*i+2]/syst->bsus_normvec[i];
		}

	}

	for (i=0;i<histogram_size-1;i++)
	{
		double lratio;
		if ( (syst->bsus_tm[3*(i+1)]>0) && (syst->bsus_tm[3*i+2]>0) )
			lratio=log(syst->bsus_tm[3*i+2]/syst->bsus_tm[3*(i+1)]);
		else
			lratio=0.;

		syst->bsus_pm[i+1]=syst->bsus_pm[i]+lratio;
	}

}


void MC_check_energy(System *syst, Output *IO) {
//void check_energy(System *syst, Output *IO) {
	int i;
	double E = 0.;
	syst->overlap = 0;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		E += MC_energy(syst, p);
		gram_schmidt(p->orientation[0], p->orientation[1], p->orientation[2]);
	}
	E *= 0.5;

	if(syst->overlap) output_exit(IO, "\nComputing energy from scratch resulted in an overlap, aborting\n");
	if(fabs(syst->energy) > 1e-5 && fabs((E - syst->energy) / syst->energy) > 1e-5) {
		output_exit(IO, "\nEnergy check failed, old energy = %lf, new energy = %lf\n", syst->energy, E);
	}

	syst->energy = E;
}
