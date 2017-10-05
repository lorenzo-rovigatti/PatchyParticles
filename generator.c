#include "defs.h"
#include "utils.h"

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

int would_overlap(System *syst, PatchyParticle *p, vector disp) {
	int ind[3], loop_ind[3];
	vector r = {p->r[0] + disp[0], p->r[1] + disp[1], p->r[2] + disp[2]};
	ind[0] = (int) ((r[0] / syst->L - floor(r[0] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[1] = (int) ((r[1] / syst->L - floor(r[1] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);
	ind[2] = (int) ((r[2] / syst->L - floor(r[2] / syst->L)) * (1. - DBL_EPSILON) * syst->cells.N_side);

	int j, k, l;
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
						vector dist = {q->r[0] - r[0], q->r[1] - r[1], q->r[2] - r[2]};
						dist[0] -= syst->L * rint(dist[0] / syst->L);
						dist[1] -= syst->L * rint(dist[1] / syst->L);
						dist[2] -= syst->L * rint(dist[2] / syst->L);

						if(SCALAR(dist, dist) < 1.) return 1;
					}
					q = q->next;
				}
			}
		}
	}

	return 0;
}

void make_initial_conf(System *syst, char *conf_name) {
	int inserted = 0;
	while(inserted < syst->N) {
		// extract new position
		vector r = { drand48() * syst->L, drand48() * syst->L, drand48() * syst->L };
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

	FILE *out = fopen(conf_name, "w");
	if(out == NULL) fprintf(stderr, "File '%s' is not writable\n", conf_name);

	fprintf(out, "0 %d %lf %lf %lf\n", syst->N, syst->L, syst->L, syst->L);

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
	new_syst.N = atoi(argv[1]);
	double density = atof(argv[2]);
	new_syst.L = pow(new_syst.N/density, 1./3.);

	Cells *cells = &new_syst.cells;
	cells->N_side = floor(new_syst.L);
	if(cells->N_side < 3) {
		cells->N_side = 3;
		fprintf(stderr, "Box side is too small, setting cells.N_side = 3\n");
	}
	cells->N = cells->N_side * cells->N_side * cells->N_side;
	cells->heads = malloc(sizeof(PatchyParticle *) * cells->N);
	fprintf(stderr, "Cells per side: %d, total: %d\n", cells->N_side, cells->N);

	int i;
	for(i = 0; i < cells->N; i++) cells->heads[i] = NULL;

	new_syst.particles = malloc(new_syst.N * sizeof(PatchyParticle));
	for(i = 0; i < new_syst.N; i++) {
		PatchyParticle *p = new_syst.particles + i;
		p->patches = NULL;
	}
	fprintf(stderr, "Creating the configuration (N = %d, L = %lf), this could take some time... ", new_syst.N, new_syst.L);

	char name[512] = "generated.rrr";
	make_initial_conf(&new_syst, name);
	fprintf(stderr, "done\n");

	free(new_syst.particles);

	return 0;
}
