#include "avb.h"
#include "MC.h"
#include "utils.h"

#include <math.h>
#include <float.h>

avbmc *avbdata;

void AVBMC_init(input_file *input, System *syst, Output *IO) {
	avbdata = malloc(sizeof(avbmc));

	int numpatches = syst->n_patches;

	avbdata->neighbours = malloc(numpatches * sizeof(PatchyParticle*));
	avbdata->num_neighbours = 0;

	avbdata->avb_vin = syst->n_patches * syst->n_patches * (M_PI * (syst->kf_delta * syst->kf_delta * syst->kf_delta + 3. * SQR(syst->kf_delta) + 3. * syst->kf_delta) * SQR(1. - syst->kf_cosmax) / 3.);

	output_log_msg(IO, "Vavb = %lf\n", avbdata->avb_vin);

	avbdata->avb_vout = syst->V - avbdata->avb_vin;
	avbdata->avb_p = 0.5;
	getInputDouble(input, "avb_p", &avbdata->avb_p, 0);

}

void AVBMC_free() {
	free(avbdata->neighbours);
	free(avbdata);
}

void _set_neighbours(System *syst, PatchyParticle *p) {
	int ind[3], loop_ind[3];
	ind[0] = (int) ((p->r[0] / syst->box[0] - floor(p->r[0] / syst->box[0])) * (1. - DBL_EPSILON) * syst->cells->N_side);
	ind[1] = (int) ((p->r[1] / syst->box[1] - floor(p->r[1] / syst->box[1])) * (1. - DBL_EPSILON) * syst->cells->N_side);
	ind[2] = (int) ((p->r[2] / syst->box[2] - floor(p->r[2] / syst->box[2])) * (1. - DBL_EPSILON) * syst->cells->N_side);

	avbdata->num_neighbours = 0;
	memset(avbdata->neighbours, 0, p->n_patches * sizeof(PatchyParticle*));

	int j, k, l, p_patch, q_patch;

	for(j = -1; j < 2; j++) {
		loop_ind[0] = (ind[0] + j + syst->cells->N_side) % syst->cells->N_side;
		for(k = -1; k < 2; k++) {
			loop_ind[1] = (ind[1] + k + syst->cells->N_side) % syst->cells->N_side;
			for(l = -1; l < 2; l++) {
				loop_ind[2] = (ind[2] + l + syst->cells->N_side) % syst->cells->N_side;
				int loop_index = (loop_ind[0] * syst->cells->N_side + loop_ind[1]) * syst->cells->N_side + loop_ind[2];

				PatchyParticle *q = syst->cells->heads[loop_index];
				while(q != NULL) {
					if(q->index != p->index) {
						int val = MC_interact(syst, p, q, &p_patch, &q_patch);

						if(val == PATCH_BOND) {
							avbdata->neighbours[p_patch] = q;
							avbdata->num_neighbours++;
						}
					}
					q = q->next;
				}
			}
		}
	}
}

void AVBMC_dynamics(System *syst, Output *IO) {
	if(drand48() < avbdata->avb_p || syst->N < 2) MC_move_rototranslate(syst, IO);
	else {
		syst->tries[AVB]++;

		PatchyParticle *receiver = syst->particles + (int) (drand48() * syst->N);

		int i;

		PatchyParticle *p;

		_set_neighbours(syst, receiver);

		// bring in
		if(drand48() > 0.5) {
			// if all the particles but rec are in rec's bonding volume then no move is possible
			if(avbdata->num_neighbours == (syst->N - 1)) return;

			int found = 0;

			// select a particle which is not in rec's bonding volume
			do {
				p = syst->particles + (int) (drand48() * syst->N);
				found = 1;
				for(i = 0; i < receiver->n_patches && found; i++)
					if(avbdata->neighbours[i] == p) found = 0;
			} while(!found || p == receiver);

			vector new_r;
			matrix new_orient;
			// target patch
			int target_patch = (int) (drand48() * p->n_patches);

			place_inside_vbonding(syst, receiver, new_r, new_orient, target_patch);

			vector disp = { new_r[0] - p->r[0], new_r[1] - p->r[1], new_r[2] - p->r[2] };

			double deltaE = -energy(syst, p);
			rototraslate_particle(syst, p, disp, new_orient);
			deltaE += energy(syst, p);

#ifdef DEBUG
			int p_patch = 0, q_patch = 0;
			assert(MC_interact(syst, receiver, p, &p_patch, &q_patch) == PATCH_BOND);
#endif

			double acc = exp(-deltaE / syst->T) * (syst->N - avbdata->num_neighbours - 1.) * avbdata->avb_vin / ((avbdata->num_neighbours + 1.) * avbdata->avb_vout);

			if(!syst->overlap && drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;
			}
			else {
				rollback_particle(syst, p);
			}
		}
		// take out
		else {
			// we need rec to have neighbors in order to take one of them out...
			if(avbdata->num_neighbours == 0) return;

			// pick a random particle from rec's neighbors
			int random_neighbour = (int) (drand48() * avbdata->num_neighbours);

			int cn = 0;
			int ci = 0;
			while(cn <= random_neighbour) {
				if(avbdata->neighbours[ci] != NULL) cn++;
				ci++;
			}

			p = avbdata->neighbours[ci - 1];

#ifdef DEBUG
			int p_patch = 0, q_patch = 0;
			assert(MC_interact(syst, receiver, p, &p_patch, &q_patch) == PATCH_BOND);
#endif

			vector new_r;
			matrix new_orient;
			vector new_patches[syst->n_patches];
			int buff;
			do {
				new_r[0] = drand48() * syst->box[0];
				new_r[1] = drand48() * syst->box[1];
				new_r[2] = drand48() * syst->box[2];
				random_orientation(syst, new_orient);
				for(i = 0; i < syst->n_patches; i++)
					MATRIX_VECTOR_MULTIPLICATION(new_orient, syst->base_patches[i], new_patches[i]);
			} while(MC_would_interact(syst, receiver, new_r, new_patches, &buff, &buff) == PATCH_BOND);

			vector disp = { new_r[0] - p->r[0], new_r[1] - p->r[1], new_r[2] - p->r[2] };

			double deltaE = -energy(syst, p);
			rototraslate_particle(syst, p, disp, new_orient);
			deltaE += energy(syst, p);

			double acc = exp(-deltaE / syst->T) * avbdata->num_neighbours * avbdata->avb_vout / ((syst->N - avbdata->num_neighbours) * avbdata->avb_vin);

			if(!syst->overlap && drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;

			}
			else {
				rollback_particle(syst, p);
			}
		}
	}

}
