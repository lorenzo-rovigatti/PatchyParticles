#include "avb.h"
#include "MC.h"
#include "output.h"
#include "parse_input.h"
#include "system.h"
#include "utils.h"

#include <math.h>
#include <float.h>

avbmc *avbdata;

void AVBMC_init(input_file *input, System *syst, Output *IO) {
	double sin_theta = sin(acos(syst->kf_cosmax));
	if(sin_theta >= 1. / (2. * (1. + syst->kf_delta))) {
		output_exit(IO, "The AVB move cannot be used with patches that do not fulfill the single-bond-per-patch condition\n");
	}

	avbdata = malloc(sizeof(avbmc));

	int numpatches = syst->n_patches;

	avbdata->neighbours = malloc(numpatches * sizeof(PatchyParticle*));
	avbdata->num_neighbours = 0;

	avbdata->avb_vin = syst->n_patches * syst->n_patches * (M_PI * (syst->kf_delta * syst->kf_delta * syst->kf_delta + 3. * SQR(syst->kf_delta) + 3. * syst->kf_delta) * SQR(1. - syst->kf_cosmax) / 3.);

	avbdata->avb_vin_per_unit=M_PI * (syst->kf_delta * syst->kf_delta * syst->kf_delta + 3. * SQR(syst->kf_delta) + 3. * syst->kf_delta) * SQR(1. - syst->kf_cosmax) / 3.;

	output_log_msg(IO, "Vavb = %lf\n", avbdata->avb_vin);

	avbdata->avb_vout = syst->V - avbdata->avb_vin;
	avbdata->avb_p = 0.5;
	/**
	 * the probability with which avb moves are attempted. The probability of regular rototranslations is just 1 - avb_p.
	 */
	getInputDouble(input, "avb_p", &avbdata->avb_p, 0);
}

void AVBMC_free() {
	free(avbdata->neighbours);
	free(avbdata);
}

void _set_neighbours(System *syst, PatchyParticle *p) {
	int ind[3], loop_ind[3];
	cells_fill_and_get_idx_from_particle(syst, p, ind);

	avbdata->num_neighbours = 0;
	memset(avbdata->neighbours, 0, p->n_patches * sizeof(PatchyParticle*));

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
							avbdata->neighbours[p_patch] = q;
							avbdata->num_neighbours++;
						}
					}
					q = syst->cells->next[q->index];
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
		_set_neighbours(syst, receiver);

		int i;
		PatchyParticle *p;

		// the move will attempt to bring p inside the bonding volume of receiver
		if(drand48() > 0.5) {
			// if all the particles but receiver are in receiver's bonding volume then no move is possible
			if(avbdata->num_neighbours == (syst->N - 1)) return;

			int found = 0;

			// select a particle which is not in receiver's bonding volume
			do {
				p = syst->particles + (int) (drand48() * syst->N);
				found = 1;
				for(i = 0; i < receiver->n_patches && found; i++)
					if(avbdata->neighbours[i] == p) found = 0;
			} while(!found || p == receiver);

			vector new_r;
			matrix new_orient;
			// choose a target patch
			int target_patch = (int) (drand48() * p->n_patches);

			place_inside_vbonding(syst, receiver, new_r, new_orient, target_patch);

			vector disp = { new_r[0] - p->r[0], new_r[1] - p->r[1], new_r[2] - p->r[2] };

			double deltaE = -syst->do_energy_old(syst, p);
			MC_rototraslate_particle(syst, p, disp, new_orient);
			deltaE += syst->do_energy_new(syst, p);

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
				MC_rollback_particle(syst, p);
			}
		}
		// the move will try to take p out out of the bonding volume of receiver
		else {
			// we need receiver to have a non-zero number of neighbors in order to take one of them out...
			if(avbdata->num_neighbours == 0) return;

			// pick a random particle from receiver's neighbors
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

			// COLORS ///////////////////////////
			// some changes here
			//vector new_r;
			matrix new_orient;
			//vector new_patches[syst->n_patches];

			int buff;
			// choose a random position and a random orientation so that, at the end of the move, p and receiver won't be bonded any more
			PatchyParticle new_particle;
			new_particle.specie=p->specie;

			do {
				new_particle.r[0] = drand48() * syst->box[0];
				new_particle.r[1] = drand48() * syst->box[1];
				new_particle.r[2] = drand48() * syst->box[2];
				random_orientation(syst, new_orient);
				for(i = 0; i < syst->n_patches; i++) {
					MATRIX_VECTOR_MULTIPLICATION(new_orient, syst->base_patches[i], new_particle.patches[i]);
				}

			//} while(MC_would_interact(syst, receiver, new_r, new_patches, &buff, &buff) == PATCH_BOND);
			} while(MC_interact(syst,receiver,&new_particle,&buff, &buff) == PATCH_BOND);


			vector disp = { new_particle.r[0] - p->r[0], new_particle.r[1] - p->r[1], new_particle.r[2] - p->r[2] };

			double deltaE = -syst->do_energy_old(syst, p);
			MC_rototraslate_particle(syst, p, disp, new_orient);
			deltaE += syst->do_energy_new(syst, p);

			double acc = exp(-deltaE / syst->T) * avbdata->num_neighbours * avbdata->avb_vout / ((syst->N - avbdata->num_neighbours) * avbdata->avb_vin);
			if(!syst->overlap && drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;

			}
			else {
				MC_rollback_particle(syst, p);
			}
		}
	}

}



void AVBMC_dynamics_colors(System *syst, Output *IO) {
	if(drand48() < avbdata->avb_p || syst->N < 2) MC_move_rototranslate(syst, IO);
	else {
		syst->tries[AVB]++;

		PatchyParticle *receiver = syst->particles + (int) (drand48() * syst->N);
		_set_neighbours(syst, receiver);

		int i;
		PatchyParticle *p;

		// the move will attempt to bring p inside the bonding volume of receiver
		if(drand48() > 0.5) {
			// if all the particles but receiver are in receiver's bonding volume then no move is possible
			if(avbdata->num_neighbours == (syst->N - 1)) return;

			int found = 0;

			// select a particle which is not in receiver's bonding volume
			do {
				p = syst->particles + (int) (drand48() * syst->N);
				found = 1;
				for(i = 0; i < receiver->n_patches && found; i++)
					if(avbdata->neighbours[i] == p) found = 0;
			} while(!found || p == receiver);

			// select which bonding volume to occupy
			int sr=receiver->specie;
			int sp=p->specie;
			int bonding_volume=syst->bonding_volume_units[sr][sp];

			// early rejection
			if (bonding_volume==0)
			{
				return;
			}

			int sel_bv=(int) (drand48() * bonding_volume);
			int p_p=-1;
			int p_r=0;
			i=-1;
			int c=syst->particlescolor[sr][p_r];
			//int cc=syst->colorint[c];

			do{

				p_p++;

				if (p_p==p->n_patches)
				{
					p_p=0;
					p_r++;

					c=syst->particlescolor[sr][p_r];
					//cc=syst->colorint[c];
				}

				int kk;
				for (kk=0;kk<syst->ncolorint[c];kk++)
				{
					int cc=syst->colorint[c][kk];

					if (syst->particlescolor[sp][p_p]==cc)
						i++;

					if (i==sel_bv)
						break;
				}

			} while (i<sel_bv);

			// early rejection
			if (avbdata->neighbours[p_r] != NULL)
			{
				return;
			}

			vector new_r;
			matrix new_orient;

			place_inside_vbonding_colored(syst, receiver, new_r, new_orient, p_r,p_p);

			vector disp = { new_r[0] - p->r[0], new_r[1] - p->r[1], new_r[2] - p->r[2] };

			double deltaE = -syst->do_energy_old(syst, p);
			MC_rototraslate_particle(syst, p, disp, new_orient);
			deltaE += syst->do_energy_new(syst, p);

#ifdef DEBUG
			int p_patch = 0, q_patch = 0;
			assert(MC_interact(syst, receiver, p, &p_patch, &q_patch) == PATCH_BOND);
#endif

			double acc = exp(-deltaE / syst->T) * (syst->N - avbdata->num_neighbours - 1.) * bonding_volume*avbdata->avb_vin_per_unit / ((avbdata->num_neighbours + 1.) * (syst->V-bonding_volume*avbdata->avb_vin_per_unit));

			if(!syst->overlap && drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;
			}
			else {
				MC_rollback_particle(syst, p);
			}
		}
		// the move will try to take p out out of the bonding volume of receiver
		else {
			// we need receiver to have a non-zero number of neighbors in order to take one of them out...
			if(avbdata->num_neighbours == 0) return;

			// pick a random particle from receiver's neighbors
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

			int sr=receiver->specie;
			int sp=p->specie;
			int bonding_volume=syst->bonding_volume_units[sr][sp];

			// COLORS ///////////////////////////
			// some changes here
			//vector new_r;
			matrix new_orient;
			vector new_patches[syst->n_patches];

			int buff;
			// choose a random position and a random orientation so that, at the end of the move, p and receiver won't be bonded any more
			PatchyParticle new_particle;
			new_particle.specie=p->specie;
			new_particle.patches=new_patches;

			do {
				new_particle.r[0] = drand48() * syst->box[0];
				new_particle.r[1] = drand48() * syst->box[1];
				new_particle.r[2] = drand48() * syst->box[2];
				random_orientation(syst, new_orient);
				for(i = 0; i < syst->n_patches; i++) {
					MATRIX_VECTOR_MULTIPLICATION(new_orient, syst->base_patches[i], new_particle.patches[i]);
				}

			//} while(MC_would_interact(syst, receiver, new_r, new_patches, &buff, &buff) == PATCH_BOND);
			} while(MC_interact(syst,receiver,&new_particle,&buff, &buff) == PATCH_BOND);


			vector disp = { new_particle.r[0] - p->r[0], new_particle.r[1] - p->r[1], new_particle.r[2] - p->r[2] };

			double deltaE = -syst->do_energy_old(syst, p);
			MC_rototraslate_particle(syst, p, disp, new_orient);
			deltaE += syst->do_energy_new(syst, p);

			double acc = exp(-deltaE / syst->T) * avbdata->num_neighbours * (syst->V-bonding_volume*avbdata->avb_vin_per_unit) / ((syst->N - avbdata->num_neighbours) * bonding_volume*avbdata->avb_vin_per_unit);
			if(!syst->overlap && drand48() < acc) {
				syst->accepted[AVB]++;
				syst->energy += deltaE;

			}
			else {
				MC_rollback_particle(syst, p);
			}
		}
	}

}
