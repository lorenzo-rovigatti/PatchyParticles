#include "vmmc.h"

#include "defs.h"
#include "MC.h"
#include "output.h"
#include "parse_input.h"
#include "system.h"
#include "utils.h"

#include <math.h>
#include <float.h>

vmmc_d *vmmcdata;

void VMMC_init(input_file *input, System *syst, Output *IO) {
	vmmcdata = malloc(sizeof(vmmc_d));

	/**
	 * if a vmmc move attempts to move a particle for more than this value, the move will be rejected
	 */
	getInputDouble(input, "vmmc_max_move", &vmmcdata->max_move, 1);
	/**
	 * if a vmmc move attempts to move more than this number of particles, the move will be rejected
	 */
	getInputInt(input, "vmmc_max_cluster", &vmmcdata->max_cluster, 1);

	vmmcdata->n_possible_links = 0;
	vmmcdata->possible_links = (PatchyParticle **) malloc(16 * syst->N_max * sizeof(PatchyParticle *)); // assuming a maximum of 8 bonds per particle

	vmmcdata->clust = (PatchyParticle **) malloc(syst->N_max * sizeof(PatchyParticle *));
	vmmcdata->n_clust = 0;

	vmmcdata->prelinked_particles = (PatchyParticle **) malloc(syst->N_max * sizeof(PatchyParticle*));
	vmmcdata->n_prelinked_particles = 0;

	vmmcdata->is_in_cluster = (int *) calloc(syst->N_max, sizeof(int*));

	vmmcdata->copies = (int *) calloc(syst->N_max, sizeof(int*));
	vmmcdata->possible_links_counter = (int *) calloc(syst->N_max, sizeof(int*));

	vmmcdata->outright_failed = (PatchyParticle **) malloc(syst->N_max * sizeof(PatchyParticle*));
	vmmcdata->n_outright_failed = 0;

	output_log_msg(IO, "Using VMMC dynamics with max_move = %lf, max_clust = %d on a system with %d particles\n", vmmcdata->max_move, vmmcdata->max_cluster, syst->N);
}

void VMMC_free() {
	free(vmmcdata->prelinked_particles);
	free(vmmcdata->clust);
	free(vmmcdata->is_in_cluster);
	free(vmmcdata->possible_links);
	free(vmmcdata->copies);
	free(vmmcdata-> possible_links_counter);
	free(vmmcdata->outright_failed);
	free(vmmcdata);
}

void _store_dof(PatchyParticle * p) {
	int i, j;
	for(i = 0; i < 3; i++) {
		p->r_old[i] = p->r[i];
		for(j = 0; j < 3; j++) {
			p->orientation_old[i][j] = p->orientation[i][j];
		}
	}
}

void _restore_dof(PatchyParticle * p) {
	int i, j, k;
	for(i = 0; i < 3; i++) {
		p->r[i] = p->r_old[i];
		for(j = 0; j < 3; j++) {
			p->orientation[i][j] = p->orientation_old[i][j];
		}
	}
	for(k = 0; k < p->n_patches; k++) {
		MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[k], p->patches[k]);
	}
}

double _pair_energy(System * syst, PatchyParticle *p, PatchyParticle *q) {
	int p_patch, q_patch;
	int val = MC_interact(syst, p, q, &p_patch, &q_patch);
	if(val == PATCH_BOND) {
		return -1.;
	}
	else if(val == OVERLAP) {
		syst->overlap = 1;
		return 0.;
	}
	else {
		return 0.;
	}
}

// TODO: maybe the following function can be avoided
int _mycomp(const void *p, const void * q, void * s) {
	PatchyParticle *a = *(PatchyParticle **) p;
	PatchyParticle *b = *((PatchyParticle **) p + 1);
	PatchyParticle *c = *(PatchyParticle **) q;
	PatchyParticle *d = *((PatchyParticle **) q + 1);
	int idx1 = a->index * ((System *) s)->N + b->index;
	int idx2 = c->index * ((System *) s)->N + d->index;
	return idx1 - idx2;
}

double _compute_cluster_energy(System *syst) {
	double res = 0.;
	int i;
	for(i = 0; i < vmmcdata->n_clust; i++) {
		PatchyParticle * p = vmmcdata->clust[i];
		assert(vmmcdata->is_in_cluster[p->index] == 1);
		int ind[3];
		cells_fill_and_get_idx_from_particle(syst, p, ind);
		int j, k, l;
		int loop_ind[3];
		for(j = -1; j < 2; j++) {
			loop_ind[0] = (ind[0] + j + syst->cells->N_side[0]) % syst->cells->N_side[0];
			for(k = -1; k < 2; k++) {
				loop_ind[1] = (ind[1] + k + syst->cells->N_side[1]) % syst->cells->N_side[1];
				for(l = -1; l < 2; l++) {
					loop_ind[2] = (ind[2] + l + syst->cells->N_side[2]) % syst->cells->N_side[2];

					int loop_index = (loop_ind[0] * syst->cells->N_side[1] + loop_ind[1]) * syst->cells->N_side[2] + loop_ind[2];
					PatchyParticle *q = syst->cells->heads[loop_index];
					while(q != NULL) {
						if(vmmcdata->is_in_cluster[q->index] == 0) {
							res += _pair_energy(syst, p, q);
						}
						q = syst->cells->next[q->index];
					}
				}
			}
		}
	}
	return res;
}

void normalized_diff_vector(System *syst, vector a, vector b, vector c) {
	c[0] = b[0] - a[0];
	c[1] = b[1] - a[1];
	c[2] = b[2] - a[2];
	c[0] -= syst->box[0] * rint(c[0] / syst->box[0]);
	c[1] -= syst->box[1] * rint(c[1] / syst->box[1]);
	c[2] -= syst->box[2] * rint(c[2] / syst->box[2]);
}

void _populate_possible_links(System *syst, Output *output_files, PatchyParticle *p, int *force_reject) {
	// get a list of possible links that can be formed by p
	int ind[3];
	cells_fill_and_get_idx_from_particle(syst, p, ind);

	assert(vmmcdata->is_in_cluster[p->index] == 1);

	int j, k, l;
	int loop_ind[3];
	for(j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + syst->cells->N_side[0]) % syst->cells->N_side[0];
		for(k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + syst->cells->N_side[1]) % syst->cells->N_side[1];
			for(l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + syst->cells->N_side[2]) % syst->cells->N_side[2];
				int loop_index = (loop_ind[0] * syst->cells->N_side[1] + loop_ind[1]) * syst->cells->N_side[2] + loop_ind[2];

				PatchyParticle *q = syst->cells->heads[loop_index];
				while(q != NULL) {
					// FIXME: following if condition may be more efficient than the one currently
					// used, but check if it creates more pain than it actually cures
					//if(vmmcdata->is_in_cluster[q->index] == 0 && q->index != p->index) {
					if(q->index != p->index) {
						vector dist = { q->r[0] - p->r[0], q->r[1] - p->r[1], q->r[2] - p->r[2] };
						dist[0] -= syst->box[0] * rint(dist[0] / syst->box[0]);
						dist[1] -= syst->box[1] * rint(dist[1] / syst->box[1]);
						dist[2] -= syst->box[2] * rint(dist[2] / syst->box[2]);
						double dist2 = SCALAR(dist, dist);

						if(dist2 < syst->kf_sqr_rcut) {
							// check that we have not exhausted the memory
							if(vmmcdata->n_possible_links >= 16 * syst->N_max) output_exit(output_files, "VMMC: memory exhausted");

							// we pick the image of q that is closest to the non virtual copy of p
							normalized_diff_vector(syst, p->r_old, q->r, dist);
							vector new_pos = { p->r_old[0] + dist[0], p->r_old[1] + dist[1], p->r_old[2] + dist[2] };

							// we check if we are using the same image of q twice for the same cluster. If we are, we early reject the move
							if(vmmcdata->copies[q->index] == 1 || vmmcdata->is_in_cluster[q->index] == 1) {
								if(fabs(new_pos[0] - q->r[0]) > 1e-6 || fabs(new_pos[1] - q->r[1]) > 1e-6 || fabs(new_pos[2] - q->r[2]) > 1e-6)	{
									*force_reject = 1;
									return;
								}
							}

							// We can update the position of q only if it has not been recruited already
							if(vmmcdata->copies[q->index] == 0 && vmmcdata->is_in_cluster[q->index] == 0) {
								q->r[0] = new_pos[0];
								q->r[1] = new_pos[1];
								q->r[2] = new_pos[2];

								vmmcdata->copies[q->index] = 1;
							}

							// add to list if entry does not exist
							vmmcdata->possible_links[2 * vmmcdata->n_possible_links + 0] = (p->index < q->index ? p : q);
							vmmcdata->possible_links[2 * vmmcdata->n_possible_links + 1] = (p->index < q->index ? q : p);
							vmmcdata->n_possible_links++;
							vmmcdata->possible_links_counter[p->index]++;
							vmmcdata->possible_links_counter[q->index]++;

							int m;
							for(m = 0; m < vmmcdata->n_possible_links - 1; m++) {
								if(_mycomp(&(vmmcdata->possible_links[2 * m]), &(vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1)]), syst) == 0) {
									vmmcdata->n_possible_links--;
									vmmcdata->possible_links_counter[p->index] --;
									vmmcdata->possible_links_counter[q->index] --;
									if(vmmcdata->possible_links_counter[p->index] == 0) vmmcdata->copies[p->index] = 0;
									if(vmmcdata->possible_links_counter[q->index] == 0) vmmcdata->copies[q->index] = 0;
									break;
								}
							}
						}
					}
					q = syst->cells->next[q->index];
				}
			}
		}
	}

	return;
}

void _move_particle(System * syst, PatchyParticle * p, vector move, double t) {
	if(vmmcdata->which_move == VMMC_TRANSLATION) {
		p->r[0] += move[0];
		p->r[1] += move[1];
		p->r[2] += move[2];
	}
	else {
		assert(vmmcdata->which_move == VMMC_ROTATION);
		// WARNING: assumes _store_dof has been called ahead of this
		// and this p->orientation_old == p_orientation
		vector dr, dr_tmp;
		PatchyParticle *seed = vmmcdata->clust[0];
		dr_tmp[0] = p->r[0] - seed->r[0];
		dr_tmp[1] = p->r[1] - seed->r[1];
		dr_tmp[2] = p->r[2] - seed->r[2];
		MATRIX_VECTOR_MULTIPLICATION(vmmcdata->rotation, dr_tmp, dr);
		p->r[0] = seed->r[0] + dr[0];
		p->r[1] = seed->r[1] + dr[1];
		p->r[2] = seed->r[2] + dr[2];
		matrix_matrix_multiplication(vmmcdata->rotation, p->orientation_old, p->orientation); // to be consistent with vector rotation, this is the correct order
		int i;
		for(i = 0; i < p->n_patches; i++) {
			MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[i], p->patches[i]);
		}
	}
}

void VMMC_dynamics(System *syst, Output *output_files) {
	syst->tries[MOVE_VMMC]++;

	// initialization of things
	int force_reject = 0;
	vmmcdata->n_clust = 1;
	vmmcdata->n_outright_failed = 0;
	vmmcdata->n_possible_links = 0;
	vmmcdata->n_prelinked_particles = 0;

	// extract particle at random and add it to cluster
	PatchyParticle *seedp = syst->particles + (int) (drand48() * syst->N);
	vmmcdata->is_in_cluster[seedp->index] = 1;
	vmmcdata->clust[0] = seedp;

	// build random move
	vector move;
	double angle = -1.;
	if(drand48() < 0.5) {
		vmmcdata->which_move = VMMC_TRANSLATION;
		move[0] = (drand48() - 0.5) * syst->disp_max;
		move[1] = (drand48() - 0.5) * syst->disp_max;
		move[2] = (drand48() - 0.5) * syst->disp_max;
	}
	else {
		vmmcdata->which_move = VMMC_ROTATION;
		random_vector_on_sphere(move);
		angle = drand48() * syst->theta_max;
		get_rotation_matrix(move, angle, vmmcdata->rotation);
	}

	// get a list of possible links before and after the move
	_store_dof(seedp);
	_populate_possible_links(syst, output_files, seedp, &force_reject);

	if(!force_reject) {
		_move_particle(syst, seedp, move, angle);
		_populate_possible_links(syst, output_files, seedp, &force_reject);
		_restore_dof(seedp);
	}

	// TODO: possibly do single particle move if seed particle has no neighbours
	while(vmmcdata->n_possible_links > 0 && vmmcdata->n_clust < vmmcdata->max_cluster && !force_reject) {
		// extract link at random from list
		int link_index = (int) (drand48() * vmmcdata->n_possible_links);

		PatchyParticle *p = vmmcdata->possible_links[2 * link_index + 0];
		PatchyParticle *q = vmmcdata->possible_links[2 * link_index + 1];

		// assert that at least one is in cluster already
		assert(vmmcdata->is_in_cluster[p->index] == 1 || vmmcdata->is_in_cluster[q->index] == 1);

		// If both are in cluster, fix the possible links array and go to process next link
		if(vmmcdata->is_in_cluster[p->index] == 1 && vmmcdata->is_in_cluster[q->index] == 1) {
			vmmcdata->possible_links[2 * link_index + 0] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 0];
			vmmcdata->possible_links[2 * link_index + 1] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 1];
			vmmcdata->n_possible_links--;
			vmmcdata->possible_links_counter[p->index] --;
			vmmcdata->possible_links_counter[q->index] --;
			if (vmmcdata->possible_links_counter[p->index] == 0) vmmcdata->copies[p->index] = 0;
			if (vmmcdata->possible_links_counter[q->index] == 0) vmmcdata->copies[q->index] = 0;
			continue;
		}

		// now we know that one particle is in the cluster and the other is not.
		// we make it so that p is in the cluster and q is not
		if(vmmcdata->is_in_cluster[p->index] == 0) {
			PatchyParticle *tmp;
			tmp = q;
			q = p;
			p = tmp;
		}

		double E_old, E_p_moved, E_q_moved;

		E_old = _pair_energy(syst, p, q);

		assert(syst->overlap == 0);

		_store_dof(p);
		_move_particle(syst, p, move, angle);
		E_p_moved = _pair_energy(syst, p, q);
		_restore_dof(p);

		int force_prelink = syst->overlap;
		syst->overlap = 0;

		double p1 = 1. - exp((1. / syst->T) * (E_old - E_p_moved));

		// decide if p wants to recruit q
		if(force_prelink == 1 || p1 > drand48()) {
			_store_dof(q);
			_move_particle(syst, q, move, angle);
			E_q_moved = _pair_energy(syst, p, q);
			_restore_dof(q);

			int force_link = syst->overlap;
			syst->overlap = 0;

			double p2 = 1. - exp((1. / syst->T) * (E_old - E_q_moved));
			if(p2 > 1.) p2 = 1.;

			// decide if q agrees to be recruited
			if(force_link == 1 || (p2 / p1) > drand48()) {
				// we recruit q
				vmmcdata->copies[q->index] = 0;
				vmmcdata->is_in_cluster[q->index] = 1;
				vmmcdata->clust[vmmcdata->n_clust] = q;
				vmmcdata->n_clust++;

				// we expand the list of possible links
				_store_dof(q);
				_populate_possible_links(syst, output_files, q, &force_reject);
				if(force_reject) break; // potential for removal, this condition is controlled by the while
				_move_particle(syst, q, move, angle);
				_populate_possible_links(syst, output_files, q, &force_reject);
				_restore_dof(q);
				if(force_reject) break; // potential for removal, this condition is controlled by the while
			}
			else {
				// if q does not want to go along with the move, it is a "prelinked particle"
				vmmcdata->prelinked_particles[vmmcdata->n_prelinked_particles] = q;
				vmmcdata->n_prelinked_particles++;
			}
		}
		else {
			vmmcdata->outright_failed[vmmcdata->n_outright_failed] = q;
			vmmcdata->n_outright_failed ++;
		}

		vmmcdata->possible_links[2 * link_index + 0] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 0];
		vmmcdata->possible_links[2 * link_index + 1] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 1];
		vmmcdata->n_possible_links--;
		vmmcdata->possible_links_counter[p->index] --;
		vmmcdata->possible_links_counter[q->index] --;
		if(vmmcdata->possible_links_counter[p->index] == 0) {
			vmmcdata->copies[p->index] = 0;
		}
		if(vmmcdata->possible_links_counter[q->index] == 0) {
			vmmcdata->copies[q->index] = 0;
		}
	}

	// we reject if the cluster is too large
	if(vmmcdata->n_clust == vmmcdata->max_cluster) force_reject = 1;

	// we reject if there are prelinked particles that have not been recruited
	int i;
	for(i = 0; i < vmmcdata->n_prelinked_particles; i++) {
		PatchyParticle * p = vmmcdata->prelinked_particles[i];
		if(vmmcdata->is_in_cluster[p->index] == 0) {
			force_reject = 1;
		}
		else {
			vmmcdata->copies[p->index] = 0;
		}
	}

	double deltaE = 0.;
	if(force_reject == 0) {
		deltaE = -_compute_cluster_energy(syst);
	}
	assert(syst->overlap == 0);

	// we move the particles and force a reject if some particle has moved too much
	// TODO: perhaps put this in the main cycle? it could be done, but putting an early rejection
	// in something as complicated as the main cycle may make things less readable.
	for(i = 0; i < vmmcdata->n_clust && !force_reject; i++) {
		PatchyParticle *p = vmmcdata->clust[i];
		_store_dof(p);
		_move_particle(syst, p, move, angle);
		MC_change_cell(syst, p);
	}

	if(force_reject == 0) {
		deltaE += _compute_cluster_energy(syst);
	}
	assert(syst->overlap == 0);

	// if we need to reject the move, we put everything back
	if(force_reject == 1) {
		for(i = 0; i < vmmcdata->n_clust; i++) {
			PatchyParticle *p = vmmcdata->clust[i];
			_restore_dof(p);
		}
	}

	// if the move is accepted, we update the simulation info
	if(force_reject == 0) {
		syst->accepted[MOVE_VMMC]++;
		syst->energy += deltaE;
	}

	// fix cells for each particle in the cluster, whether we have moved them or not
	// and we fix the values of the is_in_cluster array
	for(i = 0; i < vmmcdata->n_clust; i++) {
		PatchyParticle *p = vmmcdata->clust[i];
		MC_change_cell(syst, p);
		vmmcdata->is_in_cluster[p->index] = 0;
	}

	for(i = 0; i < vmmcdata->n_clust; i++) {
		PatchyParticle *p = vmmcdata->clust[i];
		vmmcdata->is_in_cluster[p->index] = 0;
	}

	// We reset  the copies indices of those particles that were not recruited nor prelinked
	for(i = 0; i < vmmcdata->n_outright_failed; i++) {
		PatchyParticle *p = vmmcdata->outright_failed[i];
		vmmcdata->copies[p->index] = 0;
	}
	vmmcdata->n_outright_failed = 0;


	// This for loop is entered only if we early-rejected the move because of the copies mechanism.
	// All other particles for which we have set to one the array "copies" have been reset to zero
	// as they were in the cluster or in the list of prelinked particles or in the list of outright_failed
	for(i = 0; i < vmmcdata->n_possible_links; i++) {
		PatchyParticle *p = vmmcdata->possible_links[2 * i + 0];
		PatchyParticle *q = vmmcdata->possible_links[2 * i + 1];
		vmmcdata->copies[p->index] = 0;
		vmmcdata->copies[q->index] = 0;
		vmmcdata->possible_links_counter[p->index] = 0;
		vmmcdata->possible_links_counter[q->index] = 0;
	}
}

