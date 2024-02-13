/*
 * system.c
 *
 *  Created on: 31/ott/2011
 *      Author: lorenzo
 */

#include "system.h"

#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "MC.h"
#include "output.h"
#include "utils.h"
#include "order_parameters.h"
#include "jr_cluster.h"
#include "jr_sus.h"


#define Matrix2D(array,dim1,dim2,type)                                         \
do {                                                                           \
int j;                                                                         \
array=(type**)malloc((dim1)*sizeof(type*));                                    \
array[0]=(type*)calloc((dim1)*(dim2),sizeof(type));                            \
for (j=0;j<(dim1);j++)                                                         \
array[j]=array[0]+j*(dim2);                                                    \
} while (0);



#define Free2D(array)                                                          \
free(array[0]);                                                                \
free(array);



void _init_tetrahedral_patches(System *syst, Output *output_files) {
	syst->n_patches = 4;
	syst->base_patches = malloc(sizeof(vector) * syst->n_patches);
	double half_isqrt3 = 0.5 / sqrt(3);
	set_vector(syst->base_patches[0], -half_isqrt3, -half_isqrt3,  half_isqrt3);
	set_vector(syst->base_patches[1], half_isqrt3, -half_isqrt3, -half_isqrt3);
	set_vector(syst->base_patches[2], half_isqrt3,  half_isqrt3,  half_isqrt3);
	set_vector(syst->base_patches[3], -half_isqrt3,  half_isqrt3, -half_isqrt3);

	int i, j;
	for(i = 0; i < syst->n_patches; i++) normalize(syst->base_patches[i]);

	// now we need to initialize syst->base_orient
	// first we initialize my_orient as the identity matrix
	// and we get -syst->base_patches[0], because the
	// set_orientation_around_vector invert its first argument
	matrix my_orient;
	vector my_first_patch;
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			memset(my_orient[i], 0, 3 * sizeof(double));
			my_orient[i][i] = 1.;
		}
		my_first_patch[i] = -syst->base_patches[0][i];
	}
	// then we calculate the rotation matrix required to transform
	// the 0, 0, 1 vector to the syst->base_patches[0] one
	set_orientation_around_vector(my_first_patch, my_orient, 0);
	// and then we transpose that matrix to obtain the rotation
	// needed to transform the syst->base_patches[0] vector
	// into the 0, 0, 1 one
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) syst->base_orient[i][j] = my_orient[j][i];

	}



	
}


void _init_distorted_tetrahedral_patches(System *syst, Output *output_files,double angle_a,double angle_d) {
	syst->n_patches = 4;
	syst->base_patches = malloc(sizeof(vector) * syst->n_patches);
	
	set_vector(syst->base_patches[0], -sin(angle_a)/sqrt(2), -sin(angle_a)/sqrt(2),  cos(angle_a));
	set_vector(syst->base_patches[1], sin(angle_d)/sqrt(2), -sin(angle_d)/sqrt(2), -cos(angle_d));
	set_vector(syst->base_patches[2], sin(angle_a)/sqrt(2),  sin(angle_a)/sqrt(2),  cos(angle_a));
	set_vector(syst->base_patches[3], -sin(angle_d)/sqrt(2),  sin(angle_d)/sqrt(2), -cos(angle_d));

	int i, j;
	for(i = 0; i < syst->n_patches; i++) normalize(syst->base_patches[i]);

	// now we need to initialize syst->base_orient
	// first we initialize my_orient as the identity matrix
	// and we get -syst->base_patches[0], because the
	// set_orientation_around_vector invert its first argument
	matrix my_orient;
	vector my_first_patch;
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			memset(my_orient[i], 0, 3 * sizeof(double));
			my_orient[i][i] = 1.;
		}
		my_first_patch[i] = -syst->base_patches[0][i];
	}
	// then we calculate the rotation matrix required to transform
	// the 0, 0, 1 vector to the syst->base_patches[0] one
	set_orientation_around_vector(my_first_patch, my_orient, 0);
	// and then we transpose that matrix to obtain the rotation
	// needed to transform the syst->base_patches[0] vector
	// into the 0, 0, 1 one
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) syst->base_orient[i][j] = my_orient[j][i];
	}
}


void system_init(input_file *input, System *syst, Output *output_files) {
	int res, i;

	getInputInt(input, "Dynamics", &syst->dynamics, 1);
	getInputInt(input, "Ensemble", &syst->ensemble, 1);

	// THREE_BODY ////////////
	if (getInputInt(input, "Three_body", &syst->three_body,0)==KEY_NOT_FOUND)
	{
		syst->three_body=TWO_BODY;
	}
	//////////////////////////

	getInputDouble(input, "Disp_max", &syst->disp_max, 1);
	getInputDouble(input, "Theta_max", &syst->theta_max, 1);
	getInputDouble(input, "Temperature", &syst->T, 1);

	char name[256];
	getInputString(input, "Initial_conditions_file", name, 1);

	char bsus_name[256];
	int bsus_value=getInputString(input, "Initial_bsus_file", bsus_name, 0);

	// COLORS ///////////////////////////////////////
	char colors_name[256];
	char species_name[256];
	int colors_value=getInputString(input, "Initial_colors_file", colors_name, 0);
	int species_value=getInputString(input, "Initial_species_file", species_name, 0);
	/////////////////////////////////////////////////

	if(getInputInt(input, "Seed", &syst->seed, 0) == KEY_NOT_FOUND) {
		syst->seed = time(NULL);
		output_log_msg(output_files, "Using seed %d\n", syst->seed);
	}
	srand48(syst->seed);

	FILE *conf = fopen(name, "r");
	if(conf == NULL) output_exit(output_files, "Initial_conditions_file '%s' is not readable\n", name);

	res = fscanf(conf, "%*d %d %lf %lf %lf\n", &syst->N, syst->box, syst->box + 1, syst->box + 2);
	if(res != 4) output_exit(output_files, "The initial configuration file '%s' is empty or its headers are malformed\n", name);

	getInputDouble(input, "KF_delta", &syst->kf_delta, 1);
	getInputDouble(input, "KF_cosmax", &syst->kf_cosmax, 1);

	syst->kf_sqr_rcut = SQR(1. + syst->kf_delta);

	output_log_msg(output_files, "Patch parameters: cosmax = %lf, delta = %lf\n", syst->kf_cosmax, syst->kf_delta);

	if(syst->ensemble == SUS || syst->ensemble == GC || syst->ensemble == BSUS) {
		getInputDouble(input, "Activity", &syst->z, 1);
		switch(syst->ensemble) {
		case GC:
			getInputInt(input, "GC_N_max", &syst->N_max, 1);
			syst->N_min = 0;
			break;
		case SUS:
			getInputInt(input, "Umbrella_sampling_min", &syst->N_min, 1);
			getInputInt(input, "Umbrella_sampling_max", &syst->N_max, 1);
			if(syst->N < syst->N_min) output_exit(output_files, "Number of particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->N, syst->N_min);
			if(syst->N > syst->N_max) output_exit(output_files, "Number of particles %d is larger than Umbrella_sampling_max (%d)\n", syst->N, syst->N_max);
			syst->SUS_hist = calloc(syst->N_max - syst->N_min + 1, sizeof(llint));
			break;
		case BSUS:

			getInputInt(input, "Umbrella_sampling_min", &syst->N_min, 1);
			getInputInt(input, "Umbrella_sampling_max", &syst->N_max, 1);
			if(syst->N < syst->N_min) output_exit(output_files, "Number of particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->N, syst->N_min);
			if(syst->N > syst->N_max) output_exit(output_files, "Number of particles %d is larger than Umbrella_sampling_max (%d)\n", syst->N, syst->N_max);

			int transition_size=3*(syst->N_max-syst->N_min+1);
			int histogram_size=syst->N_max-syst->N_min+1;

			syst->bsus_collect=calloc(transition_size,sizeof(double));
			syst->bsus_tm=calloc(transition_size,sizeof(double));
			syst->bsus_normvec=calloc(histogram_size,sizeof(double));
			syst->bsus_pm=calloc(histogram_size,sizeof(double));

			if (bsus_value!=KEY_NOT_FOUND)
			{
				FILE *bsus_file=fopen(bsus_name,"rb");

				if (bsus_file)
				{
					output_log_msg(output_files, "Reading initial BSUS collection matrix\n");

					char myline[512];
					int p=0;
					char *s_res = fgets(myline, 512, bsus_file);
					while(s_res != NULL) {
						sscanf(myline, "%lf %lf %lf\n", syst->bsus_collect+3*p,syst->bsus_collect+3*p+1,syst->bsus_collect+3*p+2);
						output_log_msg(output_files, "%lf %lf %lf\n",syst->bsus_collect[3*p],syst->bsus_collect[3*p+1],syst->bsus_collect[3*p+2]);

						p++;

						s_res = fgets(myline, 512, bsus_file);
					}

					assert(p==histogram_size);

					fclose(bsus_file);

					bsus_update_histo(syst);
				}
				else
				{
					output_log_msg(output_files, "No initial BSUS collection matrix found\n");
				}
			}
			else
			{
				output_log_msg(output_files, "No initial BSUS collection matrix declared\n");
			}

			break;

		default:
			output_exit(output_files, "Unsupported ensemble '%d'\n", syst->ensemble);
			break;
		}
	}
	else if (syst->ensemble == GIBBS)
	{
		getInputInt(input, "Gibbs_N_max", &syst->N_max, 1);
		syst->N_min = 0;
	}
	else {
		syst->N_max = syst->N;
	}

	syst->Lx_move = 0;
	if ((syst->ensemble == NPT) || (syst->ensemble==CNTUS)) {
		getInputDouble(input, "rescale_factor_max", &syst->rescale_factor_max, 1);
		getInputDouble(input, "P", &syst->P, 1);

		int piston_direction=-1;
		int found_pdirection=getInputInt(input,"Piston_direction",&piston_direction,0);

		if ((found_pdirection==KEY_NOT_FOUND) || (piston_direction<0) || (piston_direction>3))
		{
			syst->piston_direction=-1;
		}
		else
			syst->piston_direction=piston_direction;

	}
	else {
		getInputInt(input, "Lx_move", &syst->Lx_move, 0);
		if(syst->Lx_move) {
			//if(syst->box[1] != syst->box[2]) {
			//	output_exit(output_files, "Lx_move = 1 requires Ly = Lz\n");
			//}
			syst->yz_ratio=syst->box[1]/syst->box[2];
			
			getInputDouble(input, "Lx_change_max", &syst->Lx_change_max, 1);
			getInputDouble(input, "Lyz_min", &syst->Lyz_min, 1);
			getInputDouble(input, "Lyz_max", &syst->Lyz_max, 1);

			if(syst->Lyz_min >= syst->box[1]) {
				output_exit(output_files, "Lyz_min should be smaller than box size\n");
			}

			if(syst->Lyz_max <= syst->box[1]) {
				output_exit(output_files, "Lyz_max should be smaller than box size\n");
			}

			/* disabled relative boundaries in favor of absolute ones
			syst->Lyz_min *= syst->box[1];
			syst->Lyz_max *= syst->box[1];
			*/

			output_log_msg(output_files, "Ly and Lz are allowed to vary between %lf and %lf\n", syst->Lyz_min, syst->Lyz_max);
		}
	}


	// GIBBS initialization /////
	if(syst->ensemble == GIBBS) {
		getInputDouble(input, "Gibbs_volume_frequency", &syst->gibbsVolumeFrequency, 1);
		getInputDouble(input, "Gibbs_swap_frequency", &syst->gibbsSwapFrequency, 1);
		getInputDouble(input, "Gibbs_volume_deltamax", &syst->gibbsVolumeDeltamax, 1);
	}
	/////////////////////////////

	syst->V = syst->box[0] * syst->box[1] * syst->box[2];
	syst->particles = malloc(syst->N_max * sizeof(PatchyParticle));
	syst->energy = 0;
	syst->overlap = 0;

	int distortion;
	int distortion_value=getInputInt(input, "Distorted_tetrahedra", &distortion, 0);

	if ((distortion_value==KEY_NOT_FOUND) || (distortion==0))
	{
		output_log_msg(output_files, "Tetrahedral system\n");
		_init_tetrahedral_patches(syst, output_files);
	}
	else
	{
		double distortion_angle_donor,distortion_angle_acceptor;

		getInputDouble(input,"Distortion_angle_donor",&distortion_angle_donor,1);
		getInputDouble(input,"Distortion_angle_acceptor",&distortion_angle_acceptor,1);

		output_log_msg(output_files, "Distorted tetrahedral system with angles %lf and %lf\n",180*distortion_angle_acceptor/M_PI,180*distortion_angle_donor/M_PI);

		_init_distorted_tetrahedral_patches(syst, output_files,distortion_angle_acceptor,distortion_angle_donor);
	}

	int tot_patches=0;
	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		p->index = i;

		p->n_patches = syst->n_patches;
		p->patches = malloc(sizeof(vector) * p->n_patches);
		p->base_patches = syst->base_patches;

		if (syst->three_body==THREE_BODY)
		{
			p->patch_numneighbours=malloc(p->n_patches*sizeof(int));
			p->patch_numneighbours_old=malloc(p->n_patches*sizeof(int));
		}
		tot_patches+=p->n_patches;
	}

	/* per il potenziale a tre corpi */
	if (syst->three_body==THREE_BODY)
	{
		syst->interacting_patches=bilistaGet(tot_patches);
		syst->interacting_patches_modified=bilistaGet(tot_patches);
	}

	i = 0;
	vector p1, p2, p3;
	char myline[512];
	char *s_res = fgets(myline, 512, conf);
	while(s_res != NULL) {
		sscanf(myline, "%lf %lf %lf\n", p1, p1 + 1, p1 + 2);
		res = fscanf(conf, "%lf %lf %lf\n", p2, p2 + 1, p2 + 2);

		PatchyParticle *p = syst->particles + i;
		res = fscanf(conf, "%lf %lf %lf\n", p->r, p->r+1, p->r+2);

		// normalize the orientation matrix
		normalize(p1);
		normalize(p2);
		cross(p1, p2, p3);

		// construct the orientation matrix
		memcpy(p->orientation[0], p1, 3 * sizeof(double));
		memcpy(p->orientation[1], p2, 3 * sizeof(double));
		memcpy(p->orientation[2], p3, 3 * sizeof(double));
		gram_schmidt(p->orientation[0], p->orientation[1], p->orientation[2]);

		int j;
		for(j = 0; j < p->n_patches; j++) {
			MATRIX_VECTOR_MULTIPLICATION(p->orientation, p->base_patches[j], p->patches[j]);
		}

		i++;

		s_res = fgets(myline, 512, conf);
	}
	fclose(conf);
	if(i != syst->N) output_exit(output_files, "Number of particles found in configuration (%d) is different from the value found in the header (%d)\n", i, syst->N);


	// COLORS ////////////////////////////////////
	int num_species;
	int num_colors;
	int num_patches=syst->n_patches;

	if ((colors_value!=KEY_NOT_FOUND) && (species_value!=KEY_NOT_FOUND))
	{

		system_readColorsMax(colors_name,&num_species,&num_colors);
		syst->num_species=num_species;
		syst->num_colors=num_colors;

		Matrix2D(syst->colorint,num_colors,num_colors,int);
		syst->ncolorint=calloc(num_colors,sizeof(int));
		Matrix2D(syst->particlescolor,num_species,num_patches,int);
		Matrix2D(syst->color,num_colors,num_species,int);

		syst->species_count=calloc(num_species,sizeof(int));

		system_readColors(colors_name,syst->colorint,syst->ncolorint,syst->particlescolor,syst->color);
		system_readSpecies(species_name,syst->N,syst->num_species,syst->particles,syst->species_count);

		Matrix2D(syst->bonding_volume_units,num_species,num_species,int);

		int ii,jj;
		for (ii=0;ii<num_species;ii++)
		{
			for (jj=0;jj<num_species;jj++)
			{
				syst->bonding_volume_units[ii][jj]=0;

				int pi;
				for (pi=0;pi<syst->n_patches;pi++)
				{

					int c=syst->particlescolor[ii][pi];
					int k;
					for (k=0;k<syst->ncolorint[c];k++)
					{
						int cc=syst->colorint[c][k];
						syst->bonding_volume_units[ii][jj]+=syst->color[cc][jj];
					}
				}
			}
		}

	}
	else
	{
		num_species=1;
		num_colors=1;
		syst->num_species=1;
		syst->num_colors=1;
		Matrix2D(syst->colorint,num_colors,num_colors,int);
		syst->ncolorint=calloc(num_colors,sizeof(int));
		Matrix2D(syst->particlescolor,num_species,num_patches,int);
		Matrix2D(syst->color,num_colors,num_species,int);
		syst->species_count=calloc(num_species,sizeof(int));

		syst->colorint[0][0]=0;
		syst->ncolorint[0]=1;
		syst->particlescolor[0][0]=0;
		syst->color[0][0]=syst->n_patches;
		syst->species_count[0]=syst->N;

		int ii;
		for (ii=0;ii<syst->N;ii++)
		{
			syst->particles[ii].specie=0;
		}

		Matrix2D(syst->bonding_volume_units,num_species,num_species,int);
		syst->bonding_volume_units[0][0]+=SQR(syst->n_patches);
	}
	//////////////////////////////////////////////

	utils_reset_acceptance_counters(syst);

	syst->r_cut = 1. + syst->kf_delta;
	cells_init(syst, output_files, syst->r_cut);
	cells_fill(syst);



	// order parameters
	clustersConstructor(syst->N);
	crystalsConstructor(input,output_files,syst);

	// UMBRELLA SAMPLING ENSEMBLE
	if(syst->ensemble == CNTUS)
	{
		// NPT section
		//getInputDouble(input, "rescale_factor_max", &syst->rescale_factor_max, 1);
		//getInputDouble(input, "P", &syst->P, 1);

		// US section
		getInputInt(input,"US_bias_type",&syst->US_bias_type,1);
		getInputDouble(input,"US_harmonic_amplitude",&syst->US_k,1);
		getInputInt(input,"US_steps",&syst->US_steps,1);
		getInputDouble(input,"US_OP_0",&syst->US_OP_0,1);
		getInputDouble(input,"US_MAX_OP",&syst->US_OP_MAX,1);
		getInputDouble(input,"US_MIN_OP",&syst->US_OP_MIN,1);

		

		// allocations
		syst->US_old_pos=malloc(3*syst->N*sizeof(double));
		syst->US_old_orientation=malloc(9*syst->N*sizeof(double));

		// cluster section
		susConstructor(input,syst->N);


		int num_solid;
		syst->US_OP=(int)getOrderParameter(syst,&num_solid);

		if(syst->US_OP < syst->US_OP_MIN) output_exit(output_files, "Number of crystalline particles %d is smaller than Umbrella_sampling_min (%d)\n", syst->US_OP, syst->US_OP_MIN);
		if(syst->US_OP > syst->US_OP_MAX) output_exit(output_files, "Number of crystalline particles %d is larger than Umbrella_sampling_max (%d)\n", syst->US_OP, syst->US_OP_MAX);


		saveClusterDistribution(num_solid);

		output_log_msg(output_files,"Initial crystal size: %d\n",syst->US_OP);
		

		/*
		// copy initial conditions
		int ii;
		for (ii=0;ii<syst->N;ii++)
		{
			syst->US_old_pos[ii*3+0]=syst->particles[ii].r[0];
			syst->US_old_pos[ii*3+1]=syst->particles[ii].r[1];
			syst->US_old_pos[ii*3+2]=syst->particles[ii].r[2];

			syst->US_old_orientation[ii*9+0]=syst->particles[ii].orientation[0][0];
			syst->US_old_orientation[ii*9+1]=syst->particles[ii].orientation[0][1];
			syst->US_old_orientation[ii*9+2]=syst->particles[ii].orientation[0][2];

			syst->US_old_orientation[ii*9+3]=syst->particles[ii].orientation[1][0];
			syst->US_old_orientation[ii*9+4]=syst->particles[ii].orientation[1][1];
			syst->US_old_orientation[ii*9+5]=syst->particles[ii].orientation[1][2];

			syst->US_old_orientation[ii*9+6]=syst->particles[ii].orientation[2][0];
			syst->US_old_orientation[ii*9+7]=syst->particles[ii].orientation[2][1];
			syst->US_old_orientation[ii*9+8]=syst->particles[ii].orientation[2][2];
		}
		*/

	}


}

void system_free(System *syst) {
	cells_free(syst->cells);

	int i;
	if(syst->base_patches != NULL) free(syst->base_patches);
	for(i = 0; i < syst->N_max; i++) {
		PatchyParticle *p = syst->particles + i;
		if(p->patches != NULL) {
			free(p->patches);

			if (syst->three_body==THREE_BODY)
			{
				free(p->patch_numneighbours);
				free(p->patch_numneighbours_old);
			}
		}
	}

	if (syst->three_body==THREE_BODY)
	{
		bilistaFree(syst->interacting_patches);
		bilistaFree(syst->interacting_patches_modified);
	}

	free(syst->particles);
	if(syst->ensemble == SUS) free(syst->SUS_hist);

	if (syst->ensemble==BSUS)
	{
		free(syst->bsus_collect);
		free(syst->bsus_tm);
		free(syst->bsus_normvec);
		free(syst->bsus_pm);
	}

	if (syst->ensemble==CNTUS)
	{
		free(syst->US_old_pos);
		free(syst->US_old_orientation);
		clustersFree();
	}

	free(syst->ncolorint);
	Free2D(syst->colorint);
	Free2D(syst->particlescolor);
	Free2D(syst->color);
	Free2D(syst->bonding_volume_units);
	free(syst->species_count);

	freeCrystals();



}


void system_readColorsMax(char *namefile,int *max_species,int *max_colors)
{
	FILE *pfile=fopen(namefile,"r");

	char line[MAX_LINE_LENGTH]="";
	char buffer[MAX_LINE_LENGTH]="";

	int maxc=0;
	int maxp=0;

	while (getLine(line,pfile)>0)
	{
		if (line[0]=='B')
		{
			char *pch;
			pch = strtok (line,",");

			int c1=atoi(pch+2);

			pch = strtok (NULL,",");

			strncpy(buffer,pch,(strlen(pch)-2)*sizeof(char));
			buffer[strlen(pch)-2]='\0';

			int c2=atoi(buffer);


			if (c1>maxc)
				maxc=c1;
			if (c2>maxc)
				maxc=c2;
		}
		else if (line[0]=='C')
		{
			char *pch;
			pch = strtok (line,",");

			int p=atoi(pch+2);

			if (p>maxp)
				maxp=p;
		}

	}

	maxp++;
	maxc++;

	*max_species=maxp;
	*max_colors=maxc;

	fclose(pfile);

}


void system_readColors(char *namefile,int **colorint,int *ncolorint,int **particle,int **color)
{

	char line[MAX_LINE_LENGTH]="";
	char buffer[MAX_LINE_LENGTH]="";


	FILE *pfile=fopen(namefile,"r");

	while (getLine(line,pfile)>0)
	{
		if (line[0]=='B')
		{
			char *pch;
			pch = strtok (line,",");

			int c1=atoi(pch+2);

			pch = strtok (NULL,",");

			strncpy(buffer,pch,(strlen(pch)-2)*sizeof(char));
			buffer[strlen(pch)-2]='\0';

			int c2=atoi(buffer);

			colorint[c1][ncolorint[c1]]=c2;
			ncolorint[c1]++;

			if (c1!=c2)
			{
				colorint[c2][ncolorint[c2]]=c1;
				ncolorint[c2]++;
			}

		}
		else if (line[0]=='C')
		{
			char *pch;
			pch = strtok (line,",");

			int p=atoi(pch+2);

			pch = strtok (NULL,",");

			int pp=atoi(pch);

			pch = strtok (NULL,",");

			strncpy(buffer,pch,(strlen(pch)-2)*sizeof(char));
			buffer[strlen(pch)-2]='\0';

			int c=atoi(buffer);

			particle[p][pp]=c;

			color[c][p]++;

		}

	}

	fclose(pfile);

}


void system_readSpecies(char *nomefile,int n,int species,PatchyParticle *p,int *species_count)
{
	char line[100000]="";

	FILE *pfile=fopen(nomefile,"r");

	// intestazione
	getLine(line,pfile);
	int nn,ns;
	sscanf(line,"%d %d\n",&nn,&ns);
	//assert(n==nn);
	assert(ns==species);

	memset(species_count,0,ns*sizeof(int));

	int specie;
	int i=0;
	while (getLine(line,pfile)>0)
	{
		sscanf(line,"%d\n",&specie);

		p[i].specie=specie;
		species_count[specie]++;

		i++;
	}


	//assert(i==n);

	fclose(pfile);

}
