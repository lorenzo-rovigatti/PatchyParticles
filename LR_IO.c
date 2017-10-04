/*
 * LR_IO.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "LR_IO.h"
#include "utils.h"

void init_IO(input_file *input, LR_IO *IO) {
	IO->log = stderr;

	int restart;
	char mode[2] = "w";
	getInputInt(input, "Restart", &restart, 1);
	if(restart == 0 || restart == 1) sprintf(mode, "a");

	char name[256];
	if(getInputString(input, "Log_file", name, 0) == KEY_FOUND) {
		if(strcmp("none", name) != 0) {
			FILE *mylog = fopen(name, "w");
			if(mylog == NULL) die(IO, "Log file '%s' is not writable\n", name);
			IO->log = mylog;
		}
	}

	getInputString(input, "Energy_file", name, 1);
	IO->energy = fopen(name, mode);
	if(IO->energy == NULL) die(IO, "Energy file '%s' is not writable\n", name);

	if(getInputString(input, "Acceptance_file", name, 0) == KEY_NOT_FOUND) sprintf(name, "acc.dat");
	IO->acc = fopen(name, mode);
	if(IO->acc == NULL) die(IO, "Acceptance file '%s' is not writable\n", name);

	getInputString(input, "Configuration_prefix", IO->configuration_prefix, 1);
	getInputString(input, "Configuration_last", IO->configuration_last, 1);

	int ensemble;
	getInputInt(input, "Ensemble", &ensemble, 1);
	if(ensemble != 0) {
		getInputString(input, "Density_file", name, 1);
		IO->density = fopen(name, mode);
		if(IO->density == NULL) die(IO, "Density file '%s' is not writable\n", name);

		if(ensemble == 3 || ensemble == 6) {
			getInputString(input, "Umbrella_sampling_prefix", IO->sus_prefix, 1);
		}
	}
	else IO->density = NULL;
}

void clean_IO(LR_IO *IO) {
	fclose(IO->energy);
	fclose(IO->acc);
	if(IO->log != stderr) fclose(IO->log);
	if(IO->density != NULL) fclose(IO->density);
}

void print_output(LR_IO *IO, LR_system *syst, llint step) {
	double E = (syst->N > 0) ? syst->energy / syst->N : 0;

	fprintf(IO->energy, "%lld %lf %lf\n", step, E, syst->energy);
	fflush(IO->energy);

	if(syst->ensemble != 0) {
		fprintf(IO->density, "%lld %lf %d\n", step, syst->N / syst->V, syst->N);
		fflush(IO->density);
	}

	// acceptances
	fprintf(IO->acc, "%lld %e", step, syst->accepted[ROTO_TRASL] / (double) syst->tries[ROTO_TRASL]);
	if(syst->dynamics != 0) {
		fprintf(IO->acc, " %e", syst->accepted[AVB] / (double) syst->tries[AVB]);
	}
	if(syst->use_avb_in_in) {
		fprintf(IO->acc, " %e", syst->accepted[AVB_IN_IN] / (double) syst->tries[AVB_IN_IN]);
	}
	if(syst->ensemble != 0) {
		fprintf(IO->acc, " %e", syst->accepted[ADD]/ (double) syst->tries[ADD]);
		fprintf(IO->acc, " %e", syst->accepted[REMOVE]/ (double) syst->tries[REMOVE]);
	}

	fprintf(IO->acc, "\n");
	fflush(IO->acc);
	reset_counters(syst);

	printf("%lld %lf %lf", step, E, syst->energy);
	if(syst->ensemble != 0) printf(" %lf", syst->N / syst->V);
	printf("\n");

	_print_conf(IO, syst, step, IO->configuration_last);
}

void _print_conf(LR_IO *IO, LR_system *syst, llint step, char *name) {
	FILE *out = fopen(name, "w");
	if(out == NULL) die(IO, "File '%s' is not writable\n", name);

	fprintf(out, "%lld %lld %d 0 21\n", step, step, syst->N);
	fprintf(out, "%lf %lf %lf 0. 0. 0.2\n", syst->L, syst->L, syst->L);

	int i;
	PatchyParticle *p = syst->particles;
	for(i = 0; i < syst->N; i++) {
		fprintf(out, "%lf %lf %lf\n", p->orient[0][0], p->orient[0][1], p->orient[0][2]);
		fprintf(out, "%lf %lf %lf\n", p->orient[1][0], p->orient[1][1], p->orient[1][2]);
		fprintf(out, "%.12lf %.12lf %.12lf\n", p->r[0], p->r[1], p->r[2]);
		p++;
	}
	fclose(out);
}

void print_sus(LR_IO *IO, LR_system *syst, llint step) {
	char name[512];
	sprintf(name, "%s%lld", IO->sus_prefix, step);
	FILE *out = fopen(name, "w");
	if(out == NULL) die(IO, "SUS file '%s' is not writable\n", name);

	int i;
	for(i = 0; i < (syst->N_max - syst->N_min + 1); i++) {
		if(syst->SUS_hist[i] > 0) fprintf(out, "%d %lld\n", i + syst->N_min, syst->SUS_hist[i]);
		syst->SUS_hist[i] = 0;
	}

	fclose(out);
}

void print_e_sus(LR_IO *IO, LR_system *syst, llint step) {
	char name[512];
	sprintf(name, "%s%lld", IO->sus_prefix, step);
	FILE *out = fopen(name, "w");
	if(out == NULL) die(IO, "SUS file '%s' is not writable\n", name);

	int i, j;
	for(i = 0; i < (syst->N_max - syst->N_min + 1); i++) {
		for(j = 0; j < syst->SUS_e_bins; j++) {
			llint c = syst->SUS_e_hist[i][j];
			if(c > 0) fprintf(out, "%d %d %lld\n", i + syst->N_min, j, c);
			syst->SUS_e_hist[i][j] = 0;
		}
	}

	fclose(out);
}

void print_conf(LR_IO *IO, LR_system *syst, llint step) {
	char name[1024];
	sprintf(name, "%s%lld", IO->configuration_prefix, step);
	_print_conf(IO, syst, step, name);
	_print_conf(IO, syst, step, IO->configuration_last);
}

void log_msg(LR_IO *IO, char *format, ...) {
	va_list args;
	if(IO->log != stderr) {
		va_start(args, format);
		vfprintf(IO->log, format, args);
		va_end(args);
	}

	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	fflush(IO->log);
}

void die(LR_IO *IO, char *format, ...) {
	va_list args;
	if(IO->log != stderr) {
		va_start(args, format);
		vfprintf(IO->log, format, args);
		va_end(args);
	}

	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	fflush(IO->log);
	exit(1);
}

void die_stderr(char *format, ...) {
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);

	exit(1);
}
