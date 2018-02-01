#include "defs.h"
#include "output.h"
#include "utils.h"

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

int would_overlap(System *syst, PatchyParticle *p, vector disp) {
	int ind[3], loop_ind[3];
	vector r = {p->r[0] + disp[0], p->r[1] + disp[1], p->r[2] + disp[2]};
	cells_fill_and_get_idx_from_vector(syst, r, ind);

	int j, k, l;
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
						vector dist = {q->r[0] - r[0], q->r[1] - r[1], q->r[2] - r[2]};
						dist[0] -= syst->box[0] * rint(dist[0] / syst->box[0]);
						dist[1] -= syst->box[1] * rint(dist[1] / syst->box[1]);
						dist[2] -= syst->box[2] * rint(dist[2] / syst->box[2]);

						if(SCALAR(dist, dist) < 1.) return 1;
					}
					q = syst->cells->next[q->index];
				}
			}
		}
	}

	return 0;
}

void make_initial_conf(System *syst, char *conf_name) {
	int inserted = 0;
	while(inserted < syst->N) {
		// extract a new position
		vector r = { drand48() * syst->box[0], drand48() * syst->box[1], drand48() * syst->box[2] };
		PatchyParticle *p = syst->particles + inserted;
		p->r[0] = p->r[1] = p->r[2] = 0.;
		p->index = inserted;

		if(!would_overlap(syst, p, r)) {
			random_orientation(syst, p->orientation);
			p->r[0] = r[0];
			p->r[1] = r[1];
			p->r[2] = r[2];

			// add the particle to the new cell
			int ind[3];
			int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
			syst->cells->next[p->index] = syst->cells->heads[cell_index];
			syst->cells->heads[cell_index] = p;
			p->cell = p->cell_old = cell_index;

			inserted++;
			if(syst->N > 10 && inserted % (syst->N/10) == 0) fprintf(stderr, "Inserted %d%% of the particles (%d/%d)\n", inserted*100/syst->N, inserted, syst->N);
		}
	}

	FILE *out = fopen(conf_name, "w");
	if(out == NULL) fprintf(stderr, "File '%s' is not writable\n", conf_name);

	fprintf(out, "0 %d %lf %lf %lf\n", syst->N, syst->box[0], syst->box[1], syst->box[2]);

	int i;
	PatchyParticle *p = syst->particles;
	for(i = 0; i < syst->N; i++) {
		fprintf(out, "%lf %lf %lf\n", p->orientation[0][0], p->orientation[0][1], p->orientation[0][2]);
		fprintf(out, "%lf %lf %lf\n", p->orientation[1][0], p->orientation[1][1], p->orientation[1][2]);
		fprintf(out, "%.12lf %.12lf %.12lf\n", p->r[0], p->r[1], p->r[2]);
		p++;
	}
	fclose(out);
}

int main(int argc, char *argv[]) {
	if(argc < 3) {
		fprintf(stderr, "Usage is %s N density\n", argv[0]);
		exit(1);
	}

	System new_syst;
	new_syst.N = new_syst.N_max = atoi(argv[1]);
	double density = atof(argv[2]);
	/**
	 * It is very hard to generate very dense configurations by just randomly inserting particles. Here we
	 * set a hard maximum (rho = 0.7) above which the code will cowardly refuse to even try.
	 */
	if(density > 0.7) {
		fprintf(stderr, "It is very hard to generate very dense configurations by just randomly inserting particles. This simple generator cannot produce configurations with density higher than 0.7\n");
		exit(1);
	}
	new_syst.box[0] = pow(new_syst.N/density, 1./3.);
	new_syst.box[1] = pow(new_syst.N/density, 1./3.);
	new_syst.box[2] = pow(new_syst.N/density, 1./3.);

	Output output_files;
	output_files.log = stderr;
	cells_init(&new_syst, &output_files, 1.);

	new_syst.particles = malloc(new_syst.N * sizeof(PatchyParticle));
	int i;
	for(i = 0; i < new_syst.N; i++) {
		PatchyParticle *p = new_syst.particles + i;
		p->patches = NULL;
	}
	char name[512] = "generated.rrr";
	make_initial_conf(&new_syst, name);
	fprintf(stderr, "Generation done. The new configuration has been written to the file '%s'\n", name);

	free(new_syst.particles);

	return 0;
}
