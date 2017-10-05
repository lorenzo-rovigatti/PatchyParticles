/*
 * neighs.c
 *
 *  Created on: 08/nov/2011
 *      Author: lorenzo
 */

#include <float.h>
#include <math.h>

#include "neighs.h"
#include "output.h"
#include "utils.h"

// -1 == overlap, 0 == no interaction, 1 == isotropic bond, 2 == patch-patch bond
int kf_would_interact(System *syst, PatchyParticle *p, vector r, vector *patches, int *onp, int *onq) {
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
int kf_interact(System *syst, PatchyParticle *p, PatchyParticle *q, int *onp, int *onq) {
	return kf_would_interact(syst, p, q->r, q->patches, onp, onq);
}
