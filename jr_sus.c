#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "parse_input.h"
#include "jr_vector.h"
#include "jr_bilista.h"
#include "jr_interaction_map.h"
#include "jr_cluster.h"
#include "jr_sus.h"



#define SQR(x) ((x)*(x))
#define SCALAR(v,u) (((v)->x)*((u)->x)+((v)->y)*((u)->y)+((v)->z)*((u)->z))




// histograms
static double *Histo_Nn_biased;
static double *Histo_Nn_unbiased;
static double *Histo_op;
static double Histo_op_max;
static double Histo_op_min;
static double Histo_op_step;
static int Histo_op_size;
static int Histo_num_campionamenti;

// histograms names
static char Pn_prefix[100];
static char Nn_prefix[100];
static char Nc_unc_prefix[100];
static int US_bias_type;


void susConstructor(input_file *input,int ncolloids)
{

	// LOOKUP CONFIGURATION FILE  ////////////////////////////////////////////////////////////

	int found_pnname=getInputString(input,"Pn_prefix",Pn_prefix,0);
	int found_nnname=getInputString(input,"Nn_prefix",Nn_prefix,0);
	int found_barrname=getInputString(input,"Nc_unconstrained_prefix",Nc_unc_prefix,0);
	int found_USbias=getInputInt(input,"US_bias_type",&US_bias_type,0);


	if (found_pnname==KEY_NOT_FOUND)
		strcpy(Pn_prefix,"Pn_");

	if (found_nnname==KEY_NOT_FOUND)
		strcpy(Nn_prefix,"Nn_");

	if (found_barrname==KEY_NOT_FOUND)
		strcpy(Nc_unc_prefix,"Br_");

	if (found_USbias==KEY_NOT_FOUND)
		US_bias_type=0;

	Histo_op_min=0;
	Histo_op_max=ncolloids;
	Histo_op_step=1;



	// histograms
	Histo_Nn_biased=calloc(ncolloids,sizeof(double));
	Histo_Nn_unbiased=calloc(ncolloids,sizeof(double));
	Histo_op_size=1+(int)((Histo_op_max-Histo_op_min)/Histo_op_step);
	Histo_op=calloc(Histo_op_size,sizeof(double));
	Histo_num_campionamenti=0;
}



void freeSus()
{
	free(Histo_Nn_biased);
	free(Histo_Nn_unbiased);
	free(Histo_op);

}



void updateHistograms(double op,double op0,double k,double temperature)
{

  int Csd_size;
  double *Csd_buffer=clustersGetCsd(&Csd_size);


	int i;

	Histo_num_campionamenti++;


	int pos=(int)((op-Histo_op_min)/Histo_op_step);
#ifdef DEBUG
	assert(pos<Histo_op_size);
#endif
	Histo_op[pos]++;


	for (i=0;i<Csd_size;i++)
	{
		Histo_Nn_biased[i]+=Csd_buffer[i];

		if (US_bias_type==1)
		{
			Histo_Nn_unbiased[i]+=Csd_buffer[i]*exp(k*pow(op,2./3.)*(pow(op,1./3.)-0.5*3*pow(op0,1./3.)));
		}
		else
		{
			Histo_Nn_unbiased[i]+=Csd_buffer[i]*exp(0.5*k*SQR(op-op0)/temperature);
		}
	}
}



void saveResetHistograms(llint time)
{

  int Csd_size;
  clustersGetCsd(&Csd_size);

	if (strncmp(Pn_prefix,"none",4)!=0)
	{
		char nomefile[100];

		sprintf(nomefile,"%s%lld",Pn_prefix,time);

		FILE *pfile=fopen(nomefile,"w");

		// la prima riga indica il numero di campionamenti
		fprintf(pfile,"%d\n",Histo_num_campionamenti);

		int i;

		for (i=0;i<Histo_op_size;i++)
		{
			if (Histo_op[i]>0.5)
			{
				fprintf(pfile,"%lf %lf\n",Histo_op_min+i*Histo_op_step,Histo_op[i]);
			}
		}

		fclose(pfile);
	}

	if (strncmp(Nn_prefix,"none",4)!=0)
	{
		char nomefile[100];

		sprintf(nomefile,"%s%lld",Nn_prefix,time);

		FILE *pfile=fopen(nomefile,"w");

		// la prima riga indica il numero di campionamenti
		fprintf(pfile,"%d\n",Histo_num_campionamenti);

		int i;

		for (i=0;i<Csd_size;i++)
		{
			if (Histo_Nn_biased[i]>0.5)
			{
				fprintf(pfile,"%d %lf\n",i,Histo_Nn_biased[i]);
			}
		}

		fclose(pfile);
	}

	if (strncmp(Nc_unc_prefix,"none",4)!=0)
	{
		char nomefile[100];

		sprintf(nomefile,"%s%lld",Nc_unc_prefix,time);

		FILE *pfile=fopen(nomefile,"w");

		// la prima riga indica il numero di campionamenti
		fprintf(pfile,"%d\n",Histo_num_campionamenti);

		int i;

		for (i=0;i<Csd_size;i++)
		{
			if (Histo_Nn_unbiased[i]!=0.0)
			{
				fprintf(pfile,"%d %e\n",i,Histo_Nn_unbiased[i]);
			}
		}

		fclose(pfile);
	}

	// resettiamo gli istogrammi
	int i;

	Histo_num_campionamenti=0;
	for (i=0;i<Histo_op_size;i++)
		Histo_op[i]=0.;

	for (i=0;i<Csd_size;i++)
	{
		Histo_Nn_biased[i]=0.;
		Histo_Nn_unbiased[i]=0.;
	}


}

