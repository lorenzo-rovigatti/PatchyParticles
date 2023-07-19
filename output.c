/*
 * output.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "output.h"

#include "defs.h"
#include "MC.h"
#include "parse_input.h"
#include "utils.h"
#include "order_parameters.h"

#include <stdlib.h>
#include <unistd.h>
#include <math.h>

void output_init(input_file *input, Output *output_files) {
	output_files->log = stderr;
	output_files->restart_step_counter = 1;
	getInputInt(input, "Restart_step_counter", &output_files->restart_step_counter, 0);
	const char *mode = (output_files->restart_step_counter) ? "w" : "a";

	char name[512];
	if(getInputString(input, "Log_file", name, 0) == KEY_FOUND) {
		if(strcmp("none", name) != 0) {
			FILE *mylog = fopen(name, "a");
			if(mylog == NULL) output_exit(output_files, "Log file '%s' is not writable\n", name);
			output_files->log = mylog;
		}
	}

	/**
	 * Load the initial step from the configuration file, if the user requested to not restart the step counter
	 */
	if(output_files->restart_step_counter) {
		output_files->start_from = 0;
	}
	else {
		getInputString(input, "Initial_conditions_file", name, 1);
		FILE *conf = fopen(name, "r");
		int res = fscanf(conf, "%lld %*d %*f %*f %*f\n", &output_files->start_from);
		if(res != 1) output_exit(output_files, "Invalid initial configuration: the first value in the first row should be the time step of the configuration");
		fclose(conf);
	}

	getInputLLInt(input, "Print_every", &output_files->print_every, 1);
	getInputLLInt(input, "Save_every", &output_files->save_every, 1);

	output_files->save_also_as_mgl = 0;
	getInputInt(input, "Save_also_as_mgl", &output_files->save_also_as_mgl, 0);

	sprintf(name, "energy.dat");
	getInputString(input, "Energy_file", name, 0);
	output_files->energy = fopen(name, mode);
	if(output_files->energy == NULL) output_exit(output_files, "Energy file '%s' is not writable\n", name);

	sprintf(name, "acceptance.dat");
	getInputString(input, "Acceptance_file", name, 0);
	output_files->acc = fopen(name, mode);
	if(output_files->acc == NULL) output_exit(output_files, "Acceptance file '%s' is not writable\n", name);

	getInputString(input, "Configuration_folder", output_files->configuration_folder, 1);
	if(access(output_files->configuration_folder, W_OK) != 0) {
		output_exit(output_files, "Cannot create files in directory '%s': please make sure that the directory exists and it is writable\n", output_files->configuration_folder);
	}
	sprintf(output_files->configuration_last, "last.rrr");
	getInputString(input, "Configuration_last", output_files->configuration_last, 0);

	sprintf(output_files->specie_last, "last_specie.rrr");
	getInputString(input, "Specie_last", output_files->specie_last, 0);

	int ensemble;
	getInputInt(input, "Ensemble", &ensemble, 1);
	if(ensemble != 0) {
		sprintf(name, "density.dat");
		getInputString(input, "Density_file", name, 0);
		output_files->density = fopen(name, mode);
		if(output_files->density == NULL) output_exit(output_files, "Density file '%s' is not writable\n", name);

		if(ensemble == SUS || ensemble == BSUS) {
			getInputString(input, "Umbrella_sampling_folder", output_files->sus_folder, 1);
			if(access(output_files->sus_folder, W_OK) != 0) {
				output_exit(output_files, "Cannot create files in directory '%s': please make sure that the directory exists and it is writable\n", output_files->sus_folder);
			}
		}
	}
	else output_files->density = NULL;

	output_files->print_bonds = 0;
	getInputInt(input, "Print_bonds", &output_files->print_bonds, 0);
	if(output_files->print_bonds) {
		if(ensemble != NVT) {
			output_exit(output_files, "The printing of bonds is unavailable in non-canonical ensembles");
		}

		getInputString(input, "Bonds_folder", output_files->bonds_folder, 1);
		if(access(output_files->bonds_folder, W_OK) != 0) {
			output_exit(output_files, "Cannot create files in directory '%s': please make sure that the directory exists and it is writable\n", output_files->bonds_folder);
		}
	}

	output_files->boxshape =NULL;
	int changeshape=0;
	getInputInt(input, "Lx_move", &changeshape, 0);
	if (changeshape)
	{
		output_files->boxshape = fopen("boxshape.dat", mode);
		if(output_files->boxshape == NULL) output_exit(output_files, "boxshape.dat is not writable\n");
	}


	// order parameters
	int output_op=getInputString(input,"Order_parameter_file",name,0);

	if (output_op==KEY_FOUND)
	{
		output_files->op = fopen(name, mode);
	}
	else
		output_files->op=NULL;




}

void output_free(Output *output_files) {
	fclose(output_files->energy);
	fclose(output_files->acc);
	if(output_files->log != stderr) fclose(output_files->log);
	if(output_files->density != NULL) fclose(output_files->density);
	if(output_files->boxshape != NULL) fclose(output_files->boxshape);
}

void output_print(Output *output_files, System *syst, llint step) {
	/**
	 * Print the energy
	 */
	double E = (syst->N > 0) ? syst->energy / syst->N : 0;
	fprintf(output_files->energy, "%lld %lf\n", step, E);
	fflush(output_files->energy);

	/**
	 * Print the density, if we are simulating in a non-canonical ensemble
	 */
	if (syst->ensemble==GIBBS)
	{
		fprintf(output_files->density, "%lld %lf %d", step, syst->N / syst->V, syst->N);
		int kk;
		for (kk=0;kk<syst->num_species;kk++)
			fprintf(output_files->density, " %d",syst->species_count[kk]);
		fprintf(output_files->density, "\n");
		fflush(output_files->density);
	}
	else if(syst->ensemble != 0) {
		fprintf(output_files->density, "%lld %lf %d\n", step, syst->N / syst->V, syst->N);
		fflush(output_files->density);
	}


	/**
	 * Print the box shape
	 */
	if (syst->Lx_move)
	{
		fprintf(output_files->boxshape, "%lld %lf %lf %lf\n",step,syst->box[0],syst->box[1],syst->box[2]);
		fflush(output_files->boxshape);
	}

	/**
	 * Print acceptances for the different moves, according to the chosen dynamics
	 */
	fprintf(output_files->acc, "%lld", step);
	switch (syst->dynamics) {
	case RTMC:
		fprintf(output_files->acc, " %e", syst->accepted[ROTO_TRASL] / (double) syst->tries[ROTO_TRASL]);
		break;
	case VMMC:
		fprintf(output_files->acc, " %e", syst->accepted[MOVE_VMMC] / (double) syst->tries[MOVE_VMMC]);
		break;
	case AVBMC:
		fprintf(output_files->acc, " %e", syst->accepted[ROTO_TRASL] / (double) syst->tries[ROTO_TRASL]);
		fprintf(output_files->acc, " %e", syst->accepted[AVB] / (double) syst->tries[AVB]);
		break;
	default:
		break;
	}

	// print acceptances for ensemble moves (on the same line as above)
	switch (syst->ensemble) {
	case GC:
	case SUS:
		fprintf(output_files->acc, " %e", syst->accepted[ADD]/ (double) syst->tries[ADD]);
		fprintf(output_files->acc, " %e", syst->accepted[REMOVE]/ (double) syst->tries[REMOVE]);
		break;
	case BSUS:
		fprintf(output_files->acc, " %e", syst->accepted[ADD]/ (double) syst->tries[ADD]);
		fprintf(output_files->acc, " %e", syst->accepted[REMOVE]/ (double) syst->tries[REMOVE]);
		break;
	case NPT:
		fprintf(output_files->acc, " %e", syst->accepted[VOLUME]/ (double) syst->tries[VOLUME]);
		break;
	case GIBBS:
		fprintf(output_files->acc, " %e %e", syst->accepted[VOLUME]/ (double) syst->tries[VOLUME],syst->accepted[TRANSFER]/ (double) syst->tries[TRANSFER]);
		break;
	case CNTUS:
		fprintf(output_files->acc, " %e", syst->accepted[USCNTMOVE]/ (double) syst->tries[USCNTMOVE]);
		break;
	default:
		break;
	}

	if (syst->Lx_move) {
		fprintf(output_files->acc, " %e", syst->accepted[LX]/ (double) syst->tries[LX]);
	}


	fprintf(output_files->acc, "\n");
	fflush(output_files->acc);
	utils_reset_acceptance_counters(syst);

	/**
	 * Print energy and (if necessary) density to the standard output
	 */
	printf("%lld %lf %lf", step, E, syst->energy);
	if(syst->ensemble != 0) printf(" %lf", syst->N / syst->V);
	printf("\n");

	/**
	 * Print the current configuration
	 */
	output_save(output_files, syst, step, output_files->configuration_last);

	/**
	 * Print the bond file, if requested by the user
	 */
	if(output_files->print_bonds) {
		char name[1024];
		sprintf(name, "%s/bonds_%lld.dat", output_files->bonds_folder, step);
		output_print_bonds(output_files, syst, name);
	}


	// order parameter section
	if (output_files->op)
	{
		int num_solid;
		int size=(int)getOrderParameter(syst,&num_solid);
		fprintf(output_files->op,"%lld %d\n", step,size);
		fflush(output_files->op);
	}

}

void output_print_bonds(Output *output_files, System *syst, char *name) {
	FILE *out = fopen(name, "w");

	/**
	 * The first line of the file contains the number of particles
	 */
	fprintf(out, "%d\n", syst->N);

	int p_idx;
	int ind[3], loop_ind[3];
	for(p_idx = 0; p_idx < syst->N; p_idx++) {
		PatchyParticle *p = syst->particles + p_idx;
		int p_n_bonds = 0;
		char bond_line[512] = "";

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

							/**
							 * If p and q are bonded, we increase p's bonding counter and update the string containing the indexes
							 * of the particles p is bonded with
							 */
							if(val == PATCH_BOND) {
								p_n_bonds++;
								sprintf(bond_line, "%s%d ", bond_line, q->index);
							}
						}
						q = syst->cells->next[q->index];
					}
				}
			}
		}

		/**
		 * Print two lines per particle. The first one contains p's index and the number of its bonded neighbours
		 * The second one is just the list of p's bonded neighbours. It might be empty
		 */
		fprintf(out, "%d %d\n", p->index, p_n_bonds);
		fprintf(out, "%s\n", bond_line);
	}

	fclose(out);
}

void output_save_to_mgl(Output *output_files, System *syst, char *name) {
	FILE *out = fopen(name, "w");
	if(out == NULL) output_exit(output_files, "File '%s' is not writable\n", name);

	fprintf(out, ".Box:%lf,%lf,%lf\n", syst->box[0], syst->box[1], syst->box[2]);

	int i;
	PatchyParticle *p = syst->particles;
	for(i = 0; i < syst->N; i++) {
		vector r = {p->r[0], p->r[1], p->r[2]};
		r[0] -= syst->box[0] * floor((r[0]) / syst->box[0]);
		r[1] -= syst->box[1] * floor((r[1]) / syst->box[1]);
		r[2] -= syst->box[2] * floor((r[2]) / syst->box[2]);
		fprintf(out, "%lf %lf %lf @ 0.5 C[0.521569,0.207843,0.152941] M", r[0], r[1], r[2]);

		int j;
		for(j = 0; j < syst->n_patches; j++) {
			vector patch;
			patch[0] = (0.5+syst->kf_delta)*p->patches[j][0];
			patch[1] = (0.5+syst->kf_delta)*p->patches[j][1];
			patch[2] = (0.5+syst->kf_delta)*p->patches[j][2];
			fprintf(out, " %lf %lf %lf %lf", patch[0], patch[1], patch[2], acos(syst->kf_cosmax));
			fprintf(out, " C[0.392157,0.584314,0.929412,1]");
		}
		fprintf(out, "\n");

		p++;
	}
	fclose(out);

}

void output_save(Output *output_files, System *syst, llint step, char *name) {
	FILE *out = fopen(name, "w");
	if(out == NULL) output_exit(output_files, "File '%s' is not writable\n", name);

	fprintf(out, "%lld %d %.12lf %.12lf %.12lf\n", step, syst->N, syst->box[0], syst->box[1], syst->box[2]);

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

void output_specie_save(Output *output_files, System *syst, llint step, char *name)
{
	FILE *out = fopen(name, "w");
	if(out == NULL) output_exit(output_files, "File '%s' is not writable\n", name);

	fprintf(out,"%d %d\n",syst->N,syst->num_species);

	int i;
	for (i=0;i<syst->N;i++)
	{
		fprintf(out,"%d\n",syst->particles[i].specie);
	}
	fclose(out);

}


void output_bsus(Output *IO, System *syst, llint step) {
	char name[512];
	sprintf(name, "%s/bsus-%lld.dat", IO->sus_folder, step);
	FILE *out = fopen(name, "w");
	FILE *outlast=fopen("last_bsus.dat","w");
	FILE *outhisto=fopen("last_bsus_collect.dat","wb");
	if(out == NULL) output_exit(IO, "SUS file '%s' is not writable\n", name);

	int i;
	for(i = 0; i < (syst->N_max - syst->N_min + 1); i++) {
		// if(syst->bsus_pm[i] > 0)
		// {
			fprintf(out, "%d %lf\n", i + syst->N_min, syst->bsus_pm[i]);
			fprintf(outlast, "%d %lf\n", i + syst->N_min, syst->bsus_pm[i]);
			fprintf(outhisto,"%lf %lf %lf\n",syst->bsus_collect[3*i],syst->bsus_collect[3*i+1],syst->bsus_collect[3*i+2]);
		// }
	}

	fclose(out);
	fclose(outlast);
	fclose(outhisto);
}

/*
void output_e_sus(Output *IO, System *syst, llint step) {
	char name[512];
	sprintf(name, "%s/sus-%lld.dat", output_files->sus_folder, step);
	FILE *out = fopen(name, "w");
	if(out == NULL) output_exit(output_files, "SUS file '%s' is not writable\n", name);

	int i;
	for(i = 0; i < (syst->N_max - syst->N_min + 1); i++) {
		if(syst->SUS_hist[i] > 0) fprintf(out, "%d %lld\n", i + syst->N_min, syst->SUS_hist[i]);
		syst->SUS_hist[i] = 0;
	}

	fclose(out);
}
*/

void output_log_msg(Output *output_files, char *format, ...) {
	va_list args;
	if(output_files->log != stderr) {
		va_start(args, format);
		vfprintf(output_files->log, format, args);
		va_end(args);
	}

	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	fflush(output_files->log);
}

void output_exit(Output *output_files, char *format, ...) {
	va_list args;
	if(output_files->log != stderr) {
		va_start(args, format);
		vfprintf(output_files->log, format, args);
		va_end(args);
	}

	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	fflush(output_files->log);
	exit(1);
}

void output_exit_stderr(char *format, ...) {
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	exit(1);
}
