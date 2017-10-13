/*
 * MC.c
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#include "avb.h"
#include "vmmc.h"
#include "MC.h"
#include "utils.h"
#include "output.h"

#include <float.h>
#include <math.h>

// TODO: to be removed
void _do_widom(System *syst) {
	double widom = 0.;
	int n_widom = 0;
	PatchyParticle p;
	p.index = 100000;
	p.patches = malloc(sizeof(vector) * syst->n_patches);
	int i;
	for(i = 0; i < 1000000; i++) {
		p.r[0] = drand48() * syst->L;
		p.r[1] = drand48() * syst->L;
		p.r[2] = drand48() * syst->L;

		random_orientation(syst, p.orientation);

		int j;
		for(j = 0; j < syst->n_patches; j++) {
			MATRIX_VECTOR_MULTIPLICATION(p.orientation, syst->base_patches[j], p.patches[j]);
		}

		double E = energy(syst, &p);
		if(!syst->overlap) {
			widom += exp(-E / syst->T);
		}
		n_widom++;
	}
	widom /= n_widom;
	FILE *out = fopen("widom.dat", "a");
	fprintf(out, "%lf %lf\n", syst->N / syst->V / widom, 1. / widom);
	fclose(out);
	free(p.patches);
}

void do_NVT(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N; i++) {
		syst->do_dynamics(syst, output_files);
	}

//	if(drand48() < 0.001) _do_widom(syst);
}

void do_GC(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N_max; i++) {
		if(drand48() < 0.01) {
			MC_add_remove(syst, output_files);
		}
		else if(syst->N > 0) syst->do_dynamics(syst, output_files);
	}

//	if(drand48() < 0.001) _do_widom(syst);
}

void do_SUS(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N_max; i++) {
		if(drand48() < 0.01) {
			MC_add_remove(syst, output_files);
			syst->SUS_hist[syst->N - syst->N_min]++;
		}
		else if(syst->N > 0) syst->do_dynamics(syst, output_files);
	}
}

void MC_init(input_file *input, System *syst, Output *IO) {
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
	default:
		output_exit(IO, "Ensemble %d not supported\n", syst->ensemble);
		break;
	}

	switch(syst->dynamics) {
	case RTMC:
		syst->do_dynamics = &MC_move_rototranslate;
		break;
	case VMMC:
		vmmc_init(input, syst, IO);
		syst->do_dynamics = &vmmc_dynamics;
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
		vmmc_free();
		break;
	case AVBMC:
		AVBMC_free(syst);
		break;
	}
}

void rollback_particle(System *syst, PatchyParticle *p) {
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

	// bring back the particle in the old cell
	if(p->cell != p->cell_old) {
		if(syst->cells.heads[p->cell]->index == p->index) syst->cells.heads[p->cell] = p->next;
		else {
			PatchyParticle *q = syst->cells.heads[p->cell];
			while(q->next != p)
				q = q->next;
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

void change_cell(System *syst, PatchyParticle *p) {
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

void rototraslate_particle(System *syst, PatchyParticle *p, vector disp, vector *orient) {
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

	change_cell(syst, p);
}

// -1 == overlap, 0 == no interaction, 1 == isotropic bond, 2 == patch-patch bond
int MC_would_interact(System *syst, PatchyParticle *p, vector r, vector *patches, int *onp, int *onq) {
	vector dist = {r[0] - p->r[0], r[1] - p->r[1], r[2] - p->r[2]};
	dist[0] -= syst->L * rint(dist[0] / syst->L);
	dist[1] -= syst->L * rint(dist[1] / syst->L);
	dist[2] -= syst->L * rint(dist[2] / syst->L);

	double dist2 = SCALAR(dist, dist);

	if(dist2 < 1.) return -1;
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

// -1 == overlap, 0 == no interaction, 1 == interaction
int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, int *onp, int *onq) {
	return MC_would_interact(syst, p, q->r, q->patches, onp, onq);
}

double energy(System *syst, PatchyParticle *p) {
	syst->overlap = 0;
	double E = 0.;

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
						int val = MC_interact(syst, p, q, &p_patch, &q_patch);

						if(val == PATCH_BOND) {
							E -= 1.;
						}
						else if(val == -1) {
							syst->overlap = 1;
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
	double deltaE = -energy(syst, p);
	rototraslate_particle(syst, p, disp, new_orient);
	deltaE += energy(syst, p);

	if(!syst->overlap && (deltaE < 0. || drand48() < exp(-deltaE / syst->T))) {
		syst->energy += deltaE;
		syst->accepted[type]++;
	}
	else {
		rollback_particle(syst, p);
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

		p->r[0] = drand48() * syst->L;
		p->r[1] = drand48() * syst->L;
		p->r[2] = drand48() * syst->L;

		random_orientation(syst, p->orientation);

		int i;
		for(i = 0; i < p->n_patches; i++) {
			MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
		}

		double delta_E = energy(syst, p);
		double acc = exp(-delta_E / syst->T) * syst->z * syst->V / (syst->N + 1.);

		if(!syst->overlap && drand48() < acc) {
			syst->energy += delta_E;

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

		PatchyParticle *p = syst->particles + (int) (drand48() * syst->N);

		double delta_E = -energy(syst, p);
		double acc = exp(-delta_E / syst->T) * syst->N / (syst->V * syst->z);
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

				int i;
				for(i = 0; i < 3; i++) {
					memcpy(p->orientation[i], q->orientation[i], sizeof(double) * 3);
				}
				for(i = 0; i < p->n_patches; i++) {
					MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
				}

				// and finally add it back to its former cell
				p->cell = q->cell;
				p->next = syst->cells.heads[p->cell];
				syst->cells.heads[p->cell] = p;
			}

			syst->accepted[REMOVE]++;
		}
	}
}

void check_energy(System *syst, Output *IO) {
	int i;
	double E = 0.;
	syst->overlap = 0;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		E += energy(syst, p);
		gram_schmidt(p->orientation[0], p->orientation[1], p->orientation[2]);
	}
	E *= 0.5;

	if(syst->overlap) output_exit(IO, "\nComputing energy from scratch resulted in an overlap, aborting");
	if(fabs(syst->energy) > 1e-5 && fabs((E - syst->energy) / syst->energy) > 1e-5) {
		output_exit(IO, "\nEnergy check failed, old energy = %lf, new energy = %lf\n", syst->energy, E);
	}

	syst->energy = E;
}
