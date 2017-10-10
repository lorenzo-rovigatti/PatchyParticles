#include "vmmc.h"
#include "MC.h"
#include "utils.h"
#include "neighs.h"

#include <math.h>
#include <float.h>

vmmc_d * vmmcdata;

void vmmc_init(input_file *input, System *syst, Output *IO) {

	vmmcdata = malloc(sizeof(vmmc_d));
	
	getInputDouble(input, "vmmc_max_move", &vmmcdata->max_move, 1);
	// TODO: fix vmmc_max_cluster to N and make option non-mandatory
	getInputInt(input, "vmmc_max_cluster", &vmmcdata->max_cluster, 1);
	
	vmmcdata->n_possible_links = 0;
	vmmcdata->possible_links = (PatchyParticle **) malloc (16 * syst->N * sizeof(PatchyParticle *)); // assuming a maximum of 8 bonds per particle
	// TODO: implement a fix if we get too many links

	vmmcdata->clust = (PatchyParticle **) malloc(syst->N * sizeof(PatchyParticle *));
	vmmcdata->n_clust = 0;
	
	vmmcdata->prelinked_particles = (PatchyParticle **) malloc(syst->N * sizeof(PatchyParticle*));
	vmmcdata->n_prelinked_particles = 0;
	
	vmmcdata->is_in_cluster = (int *) calloc(syst->N, sizeof(int*));

	output_log_msg(IO, "Using VMMC dynamics with max_move = %lf, max_clust = %d on a system with %d particles\n", vmmcdata->max_move, vmmcdata->max_cluster, syst->N);
}

void vmmc_free() {
	free(vmmcdata->prelinked_particles);
	free(vmmcdata->clust);
	free(vmmcdata->is_in_cluster);
	free(vmmcdata->possible_links);
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
	int i, j;
	for(i = 0; i < 3; i++) {
		p->r[i] = p->r_old[i];
		for(j = 0; j < 3; j++) {
			p->orientation[i][j] = p->orientation_old[i][j];
		}
	}
}

double _pair_energy(System * syst, PatchyParticle *p, PatchyParticle *q) {
	int p_patch, q_patch;
	int val = kf_interact(syst, p, q, &p_patch, &q_patch);
	if(val == PATCH_BOND) {
		return -1.;
	}
	else if(val == -1) {
		syst->overlap = 1;
		return 0.;
	}
	else {
		return 0.;
	}
}

int _mycomp (const void *p, const void * q, void * s) {
	PatchyParticle * a = (PatchyParticle *) p; 
	PatchyParticle * b = (PatchyParticle *) p + 1; 
	PatchyParticle * c = (PatchyParticle *) q; 
	PatchyParticle * d = (PatchyParticle *) q + 1; 
	int idx1 = a->index * ((System *)s)->N + b->index;
	int idx2 = c->index * ((System *)s)->N + d->index;
	printf ("N, idx1, idx2: %d, %d, %d\n", ((System *)s)->N, idx1, idx2);
	if (idx1 == idx2) return -1;
	else return idx1 - idx2;
}

void _populate_possible_links(System * syst, PatchyParticle *p) {
	// get a list of possible links that can be formed by p
	int ind[3];
	ind[0] = (int) ((p->r[0] / syst->L - floor(p->r[0] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[1] = (int) ((p->r[1] / syst->L - floor(p->r[1] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[2] = (int) ((p->r[2] / syst->L - floor(p->r[2] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	
	int j, k, l;
	int loop_ind[3];
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
						vector dist = {q->r[0] - p->r[0], q->r[1] - p->r[1], q->r[2] - p->r[2]};
						dist[0] -= syst->L * rint(dist[0] / syst->L);
						dist[1] -= syst->L * rint(dist[1] / syst->L);
						dist[2] -= syst->L * rint(dist[2] / syst->L);

						double dist2 = SCALAR(dist, dist);
						if (dist2 < syst->kf_sqr_rcut) {
							// add to list
							vmmcdata->possible_links[2 * vmmcdata->n_possible_links + 0] = (p->index < q->index ? p : q);
							vmmcdata->possible_links[2 * vmmcdata->n_possible_links + 1] = (p->index < q->index ? q : p);
							vmmcdata->n_possible_links ++;
						}
					}
					q = q->next;
				}
			}
		}
	}
	
	int i;
	printf ("@@BEGIN ## \n" );
	for (i= 0; i < vmmcdata->n_possible_links; i ++) {
		PatchyParticle * u = vmmcdata->possible_links[2 * i + 0];
		PatchyParticle * v = vmmcdata->possible_links[2 * i + 1];
		printf ("@@## %d %d\n", u->index, v->index);
	}
	printf ("@@ NOW SORT## \n" );
	
	qsort_r (vmmcdata->possible_links, vmmcdata->n_possible_links, 2 * sizeof(PatchyParticle *), _mycomp, (void *)syst);

	for (i= 0; i < vmmcdata->n_possible_links; i ++) {
		PatchyParticle * u = vmmcdata->possible_links[2 * i + 0];
		PatchyParticle * v = vmmcdata->possible_links[2 * i + 1];
		printf ("@@ ## %d %d\n", u->index, v->index);
	}
	printf ("@@ SORTED ## \n" );
	
	// remove those horrible duplicates
	int check = 1;
	while (check != 0) {
		check = 0;
		int idx1, idx2;
		for (i = 0; i < vmmcdata->n_possible_links && check == 0; i ++) {
			PatchyParticle * u = vmmcdata->possible_links[2 * i + 0];
			PatchyParticle * v = vmmcdata->possible_links[2 * i + 1];
			idx1 = u->index * syst->N + v->index;
			int j;
			for (j = 0; j < i && check == 0; j ++) {
				PatchyParticle * x = vmmcdata->possible_links[2 * j + 0];
				PatchyParticle * y = vmmcdata->possible_links[2 * j + 1];
				idx2 = x->index * syst->N + y->index;
				if (idx1 == idx2) {
					vmmcdata->possible_links[2*i
					// FIXME!!!!
			}
		}
	}
	return;
}

void _move_particle(PatchyParticle * p, vector move, double t) {
	if (vmmcdata->which_move == VMMC_TRANSLATION) {
		p->r[0] += move[0];
		p->r[1] += move[1];
		p->r[2] += move[2];
	}
	else {
		// WARNING: assumes store_dof has been called correctly
		utils_rotate_matrix(p->orientation_old, p->orientation, move, t);
	}
}

void vmmc_dynamics(System *syst, Output *IO) {
	syst->tries[MOVE_VMMC] ++;
	
	// initialization of things
	vmmcdata->n_possible_links = 0;
	vmmcdata->n_clust = 1;
	vmmcdata->n_prelinked_particles = 0;
	
	// extract particle at random
	PatchyParticle *seedp = syst->particles + (int) (drand48() * syst->N);
	vmmcdata->is_in_cluster[seedp->index] = 1;
	
	vector move;
	double angle = -1.;
	
	// random move
	if (drand48() < 0.5) {
		vmmcdata->which_move = VMMC_TRANSLATION;
		move[0] = (drand48() - 0.5) * syst->disp_max;
		move[1] = (drand48() - 0.5) * syst->disp_max;
		move[2] = (drand48() - 0.5) * syst->disp_max;
	}
	else {
		vmmcdata->which_move = VMMC_ROTATION;
		random_vector_on_sphere(move);
		angle = drand48() * syst->theta_max;
	}
	
	// get a list of possible links
	_populate_possible_links(syst, seedp);
	
	_store_dof(seedp);
	_move_particle(seedp, move, angle);
	_populate_possible_links(syst, seedp);
	_restore_dof(seedp);

	// TODO: possibly do single particle move if seed particle has no neighbours
	
	while (vmmcdata->n_possible_links > 0 && vmmcdata->n_clust < vmmcdata->max_cluster) {
		// extract link at random
		int link_index = (int) (drand48() * vmmcdata->n_possible_links);
		
		PatchyParticle * p = vmmcdata->possible_links[2 * link_index + 0];
		PatchyParticle * q = vmmcdata->possible_links[2 * link_index + 1];
		
		printf("## p->index, q->index, %d, %d\n", p->index, q->index);
		
		// assert that at least one is in cluster already
		assert (vmmcdata->is_in_cluster[p->index] == 1 || vmmcdata->is_in_cluster[q->index] == 1);
		
		// if both are in cluster, fix the possible links array and exit now
		if (vmmcdata->is_in_cluster[p->index] == 1 && vmmcdata->is_in_cluster[q->index] == 1) {
			vmmcdata->possible_links[2 * link_index + 0] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 0];
			vmmcdata->possible_links[2 * link_index + 1] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 1];
			vmmcdata->n_possible_links --;
			continue;
		}
		
		// now we know that one particle is in the cluster and the other is not.
		// we make it so that p is in the cluster and q is not
		if (vmmcdata->is_in_cluster[p->index] != 1) {
			PatchyParticle *tmp;
			tmp = q;
			q = p;
			p = tmp;
		}
	
		double E_old, E_p_moved, E_q_moved;
		
		E_old = _pair_energy(syst, p, q);
		assert (syst->overlap == 0);
		
		_store_dof(p);
		_move_particle(p, move, angle);
		E_p_moved = _pair_energy(syst, p, q);
		_restore_dof(p);
		
		int force_prelink = syst->overlap;
		syst->overlap = 0;
		
		// TODO: remove this comment; number VMMC_link(double E_new, double E_old) { return (1. - exp((1. / this->_T) * (E_old - E_new)));}
		double p1 = 1. - exp((1. / syst->T) * (E_old - E_p_moved));
		
		if (force_prelink == 1 || p1 > drand48()) {
			_store_dof(q);
			_move_particle(q, move, angle);
			E_q_moved = _pair_energy(syst, p, q);
			_restore_dof(q);
			int force_link = syst->overlap;
			syst->overlap = 0;
			double p2 = 1. - exp((1. / syst->T) * (E_old - E_q_moved));
			if (p2 > 1.) p2 = 1.;

			if (force_link == 1 || (p2 / p1) > drand48()) {
				// we recruit q
				vmmcdata->is_in_cluster[q->index] = 1;
				vmmcdata->clust[vmmcdata->n_clust] = q;
				vmmcdata->n_clust ++;
				
				// we expand the list of possible links
				_populate_possible_links(syst, q);
				
				_store_dof(q);
				_move_particle(q, move, angle);
				_populate_possible_links(syst, q);
				_restore_dof(q);
			}
			else {
				vmmcdata->prelinked_particles[vmmcdata->n_prelinked_particles] = q;
				vmmcdata->n_prelinked_particles ++;
			}
		}
		
		vmmcdata->possible_links[2 * link_index + 0] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 0];
		vmmcdata->possible_links[2 * link_index + 1] = vmmcdata->possible_links[2 * (vmmcdata->n_possible_links - 1) + 1];
		vmmcdata->n_possible_links --;
	}

	int force_reject = 0;
	// we reject if the cluster is too large
	if (vmmcdata->n_clust == vmmcdata->max_cluster)
		force_reject = 1;
	// we reject if there are prelinked particles that are not in the cluster
	int i;
	for (i = vmmcdata->n_prelinked_particles - 1; i >=0 && force_reject == 0; i --) {
		PatchyParticle * p = vmmcdata->prelinked_particles[i];
		if (vmmcdata->is_in_cluster[p->index] == 0)
			force_reject = 1;
	}
	// we move the particles and force a reject if sombody moved too much
	// TODO: put this in the main cycle
	for (i = 0; i < vmmcdata->n_clust && force_reject == 0; i ++) {
		PatchyParticle *p = vmmcdata->clust[i];
		vector old_pos = {p->r[0], p->r[1], p->r[2]};
		_store_dof(p);
		_move_particle(p, move, angle);
		vector dist = {old_pos[0] - p->r[0], old_pos[1] - p->r[1], old_pos[2] - p->r[2]};
		double dist2 = SCALAR(dist, dist);
		if (dist2 > vmmcdata->max_move * vmmcdata->max_move) {
			force_reject = 1;
		}
	}
	
	if (force_reject == 1) {
		for (i = 0; i < vmmcdata->n_clust && force_reject == 0; i ++) {
			PatchyParticle *p = vmmcdata->clust[i];
			_restore_dof(p);
		}
	}
	
	if (force_reject == 0)
		syst->accepted[MOVE_VMMC] ++;

	// fix cells for each particle in the cluster, whether we have moved them or not
	for (i = 0; i < vmmcdata->n_clust; i ++) {
		PatchyParticle *p = vmmcdata->clust[i];
		change_cell(syst, p);
		vmmcdata->is_in_cluster[p->index] = 0;
	}
	
	// checks
	// TODO: remove
	int check = 0;
	for (i = 0; i < syst->N; i ++) {
		check += vmmcdata->is_in_cluster[i];
	}
	assert (check = 0);
	
	return;
}


