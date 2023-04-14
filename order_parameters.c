#include "order_parameters.h"
#include "MC.h"
#include "output.h"
#include "parse_input.h"
#include "system.h"
#include "utils.h"





#include "jr_interaction_map.h"
#include "jr_cell_list.h"
#include "jr_smart_allocator.h"
#include "jr_cluster.h"
#include "jr_bilista.h"
#include "jr_ordinator.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>





static double Range;
static int MaxNeighbours;
static int Method;

// for Method 0
static double Coherence_threshold;
static int Solid_threshold;
// for Method 3
static double Staggered_threshold;
static double Eclipsed_lowthreshold;
static double Eclipsed_highthreshold;

static int OP_l;
static char OP_code;
static orderparam *OP;
static complex double ***Punt_qlm;


static interactionmap *Coherent_map;
static int *Num_coherent;

static double *Buffer;


static bilista *List_solid;
static int Num_solid=0;
static clusters *Cluster_solid;

static double *Csd_buffer;
static int Csd_size;

// reference system
static vector Z_versor;
static vector X_versor;
static vector Y_versor;



static listcell *CrystalFractionCell;
static int *ImSolid;


static sph_ws Legendre_workspace;
static double *Legendre_associated_polynomials;

// FUNZIONI PER IL CALCOLO VELOCE DELLE ARMONICHE SFERICHE

static double intpow (const double x, const int i) {
	if (i < 0) abort();
	if (i == 0) return 1.;
	if (i == 1) return x;
	else return x * intpow(x, i - 1);
}

static double cintpow (const complex double x, const int i) {
	if (i < 0) abort();
	if (i == 0) return 1.;
	if (i == 1) return x;
	else return x * cintpow(x, i - 1);
}


static double binomial (const double n, const int k) {
	int i;
   	double res = 1.;
	for (i = 1; i <= k; i ++) res *= (n - (k - i)) / (double) i;
	return res;
}



static double coeff (int l, int m) {
	int r = (l - m) / 2;
	int d = (l - m) % 2;
	if (d != 0) return 0.;
	else return (r%2?-1.:1.) * intpow (1./2., l) * binomial (2 * l - 2 * r, l) * binomial (l, r);
}

unsigned long long int fact(unsigned long long int n)
{
	if (n<=1) return 1;
	else return n*fact(n-1);
}

void sph_ws_init (sph_ws * w, const int lmax) {
	int l, m;

	w->lmax = lmax;

	w->PL0 = (double**) calloc ((lmax + 1), sizeof (double * ));
	w->nrm = (double**) calloc ((lmax + 1), sizeof (double * ));
	for (l = 0; l <= lmax; l ++) {
		w->PL0[l] = (double *) calloc ((l + 1), sizeof (double));
		w->nrm[l] = (double *) calloc ((l + 1), sizeof (double));
	}
	w->tmp = (double *) calloc ((lmax + 2), sizeof (double)); // we need extra space later
	w->xp = calloc ((l + 1), sizeof (double));

	// now we set the initial things
	// http://en.wikipedia.org/wiki/Legendre_polynomials
	for (l = 0; l <= lmax; l ++) {
		for (m = 0; m <= l; m ++) {
			w->PL0[l][m] = coeff (l, m);
			w->nrm[l][m] = sqrt ((2*l+1.)/ (4.*M_PI)) * sqrt (tgamma(l-m+1.)) / sqrt(tgamma(l+m+1.)) ;
		}
	}
}

void sph_ws_free (sph_ws * w) {
	int i;
	for (i = 0; i <= w->lmax; i ++) {
		free (w->PL0[i]);
		free (w->nrm[i]);
	}
	free (w->PL0);
	free (w->tmp);
	free (w->xp);
}



void dPl0dm (sph_ws * w, double * res, const double x, const int l) {
	int k, m;

	w->xp[0] = 1;
	for (m = 1; m <=l; m ++) {
		w->xp[m] = x * w->xp[m-1];
	}

	for (k = 0; k <= l; k ++) w->tmp[k] = w->PL0[l][k];

	double seno=sqrt(1. - x*x);
	double seno_last=1;

	m = 0;
	for (m = 0; m <= l; m ++) {
		res[m] = 0.;
		for (k = 0; k <= l - m; k ++) {
			res[m] += w->tmp[k] * w->xp[k];
			w->tmp[k] = (k + 1) * w->tmp[k + 1];
		}
		//res[m] *= intpow(-1, m) * pow(1. - x*x, m / 2.);
		res[m] *= (m % 2 ? -1. : 1.) * seno_last*w->nrm[l][m];
		seno_last*=seno;
	}
}

orderparam *allocateOP(int ncolloids,int MAX_LS)
{
	orderparam *op=malloc(sizeof(orderparam));

	op->o=getOrdinator(MAX_LS);

	//op->im=createInteractionMap(ncolloids,ncolloids);
	op->ime=createInteractionMap(ncolloids,ncolloids);

	op->q=calloc(MAX_LS,sizeof(complex double**));
	op->q_list=calloc(MAX_LS,sizeof(int));
	op->q_num=0;

	op->Q=calloc(MAX_LS,sizeof(complex double**));
	op->Q_list=calloc(MAX_LS,sizeof(int));
	op->Q_num=0;


	op->ql=calloc(MAX_LS,sizeof(double*));

	op->Ql=calloc(MAX_LS,sizeof(double*));

	op->wl=calloc(MAX_LS,sizeof(double*));

	op->Wl=calloc(MAX_LS,sizeof(double*));


	return op;
}

void resetOP(orderparam *op,int ncolloids,int MAX_LS)
{
	// questo reset non prevede la possibilita' di fare altre request
	// serve per quando nel corso di una simulazione si vogliono
	// calcolare tante volte le stesse quantita'

	int l,m,i;

	//resetInteractionMap(op->im);
	resetInteractionMap(op->ime);


	int il,order;
	for (il=0;il<op->q_num;il++)
	{
		l=op->q_list[il];

		order=ordinatorGetOrder(op->o,l);

		for (i=0;i<ncolloids;i++)
				for (m=-l;m<=l;m++)
					op->q[order][i][m+l]=0.;
	}

	for (il=0;il<op->Q_num;il++)
	{
		l=op->Q_list[il];

		order=ordinatorGetOrder(op->o,l);

		for (i=0;i<ncolloids;i++)
				for (m=-l;m<=l;m++)
					op->Q[order][i][m+l]=0.;
	}

	for (l=0;l<MAX_LS;l++)
	{

		/*
		 * non necessari in quanto non hanno bisogno di strutture dati inizializzate
		if (op->ql[l]!=NULL)
		{
			for (i=0;i<ncolloids;i++)
				op->ql[l][i]=0.;
		}

		if (op->Ql[l]!=NULL)
		{
			for (i=0;i<ncolloids;i++)
				op->Ql[l][i]=0.;
		}
		*/

		if (op->wl[l]!=NULL)
		{
			for (i=0;i<ncolloids;i++)
				op->wl[l][i]=0.;
		}

		if (op->Wl[l]!=NULL)
		{
			for (i=0;i<ncolloids;i++)
				op->Wl[l][i]=0.;
		}

	}

}

void freeOP(orderparam *op,int MAX_LS)
{
	int i;

	//freeInteractionMap(op->im);
	freeInteractionMap(op->ime);
	freeOrdinator(op->o);

	for (i=0;i<MAX_LS;i++)
	{
		if (op->q[i]!=NULL)
			Free2D(op->q[i]);
		if (op->Q[i]!=NULL)
			Free2D(op->Q[i]);
		if (op->ql[i]!=NULL)
			free(op->ql[i]);
		if (op->Ql[i]!=NULL)
			free(op->Ql[i]);
		if (op->wl[i]!=NULL)
			free(op->wl[i]);
		if (op->Wl[i]!=NULL)
			free(op->Wl[i]);

	}

	free(op->q);
	free(op->q_list);

	free(op->Q);
	free(op->Q_list);

	free(op->ql);

	free(op->Ql);

	free(op->wl);

	free(op->Wl);

	free(op);
}

int request_qlm(int l,int ncolloids,orderparam *op)
{
	int order=ordinatorGetOrder(op->o,l);

	if (op->q[order]==NULL)
	{
		complex double **q;
		Matrix2D(q,ncolloids,2*l+1,complex double);
		op->q[order]=q;

		op->q_list[op->q_num++]=l;
	}

	return order;
}

int request_Qlm(int l,int ncolloids,orderparam *op)
{
	// dependency
	int order=request_qlm(l,ncolloids,op);

	if (op->Q[order]==NULL)
	{
		complex double **Q;
		Matrix2D(Q,ncolloids,2*l+1,complex double);
		op->Q[order]=Q;

		op->Q_list[op->Q_num++]=l;
	}

	return order;
}

int request_ql(int l,int ncolloids,orderparam *op)
{
	// dependency
	int order=request_qlm(l,ncolloids,op);

	if (op->ql[order]==NULL)
	{
		double *ql=calloc(ncolloids,sizeof(double));
		op->ql[order]=ql;
	}

	return order;
}

int request_Ql(int l,int ncolloids,orderparam *op)
{
	// dependency
	int order=request_Qlm(l,ncolloids,op);

	if (op->Ql[order]==NULL)
	{
		double *Ql=calloc(ncolloids,sizeof(double));
		op->Ql[order]=Ql;
	}

	return order;
}

int request_wl(int l,int ncolloids,orderparam *op)
{
	// dependency
	int order=request_qlm(l,ncolloids,op);

	if (op->wl[order]==NULL)
	{
		double *wl=calloc(ncolloids,sizeof(double));
		op->wl[order]=wl;
	}

	return order;
}

int request_Wl(int l,int ncolloids,orderparam *op)
{
	// dependency
	int order=request_Qlm(l,ncolloids,op);

	if (op->Wl[order]==NULL)
	{
		double *Wl=calloc(ncolloids,sizeof(double));
		op->Wl[order]=Wl;
	}

	return order;
}


void calculate_qlm_maxneighbours(orderparam *op,vector *pos,vector *Box,double range,int maxneighbours,listcell *cells,vector x_versor,vector y_versor,vector z_versor)
{
	interactionmap *ime=op->ime;

	calculateInteractionMapWithCutoffDistanceOrdered(cells,ime,pos,Box,range);

	// settiamo il massimo numero di neighbours
	int index;
	for (index=0;index<ime->num;index++)
	{
		//printf("%d\n",ime->howmany[index]);
		if (ime->howmany[index]>maxneighbours)
			ime->howmany[index]=maxneighbours;
	}

	vector v_12,xy_projection,z_projection;


	// scegliamo una terna di riferimento casuale in modo da poter partire anche dal cristallo

	//vector old_v_12;
	double inorm;
	double costeta,cosphi,sinphi,phi;
	int particle1,particle2,j,l,m;

	for (particle1=0;particle1<ime->num;particle1++)
	{

		for (j=0;j<ime->howmany[particle1];j++)
		{
			particle2=ime->with[particle1][j];

			// calcoliamo il versore

			v_12.x=pos[particle1].x-pos[particle2].x;
			v_12.y=pos[particle1].y-pos[particle2].y;
			v_12.z=pos[particle1].z-pos[particle2].z;

			v_12.x-=Box->x*rint(v_12.x/Box->x);
			v_12.y-=Box->y*rint(v_12.y/Box->y);
			v_12.z-=Box->z*rint(v_12.z/Box->z);


			inorm=1./sqrt(SQR(v_12.x)+SQR(v_12.y)+SQR(v_12.z));

			v_12.x*=inorm;
			v_12.y*=inorm;
			v_12.z*=inorm;

			costeta=SCALAR(&v_12,&z_versor);


			double z_projection_module=SCALAR(&v_12,&z_versor);

			z_projection.x=z_projection_module*z_versor.x;
			z_projection.y=z_projection_module*z_versor.y;
			z_projection.z=z_projection_module*z_versor.z;

			xy_projection.x=v_12.x-z_projection.x;
			xy_projection.y=v_12.y-z_projection.y;
			xy_projection.z=v_12.z-z_projection.z;

			inorm=1./sqrt(SQR(xy_projection.x)+SQR(xy_projection.y)+SQR(xy_projection.z));
			xy_projection.x*=inorm;
			xy_projection.y*=inorm;
			xy_projection.z*=inorm;

			cosphi=SCALAR(&xy_projection,&x_versor);
			sinphi=SCALAR(&xy_projection,&y_versor);

			if (sinphi>0)
				phi=acos(cosphi);
			else
				phi=2.*M_PI-acos(cosphi);

			int il;
			for (il=0;il<op->q_num;il++)
			{
				l=op->q_list[il];

				int order=ordinatorGetOrder(op->o,l);

				for (m=-l;m<=l;m++)
				{

					complex double spherical_harmonic;
					if (m>=0)
					{
						spherical_harmonic=gsl_sf_legendre_sphPlm(l,m,costeta)*cexp(I*m*phi);
					}
					else
					{
						spherical_harmonic=pow(-1.,-m)*gsl_sf_legendre_sphPlm(l,-m,costeta)*cexp(I*m*phi);
					}


					op->q[order][particle1][m+l]+=spherical_harmonic/(ime->howmany[particle1]);

					// ATTENZIONE MODIFICA IMPORTANTE
					//op->q[order][particle2][m+l]+=pow(-1,l)*spherical_harmonic/(ime->howmany[particle2]);
				}
			}

		}

	}

}


void calculate_qlm_maxneighbours_fast(orderparam *op,vector *pos,vector *Box,double range,int maxneighbours,listcell *cells,vector x_versor,vector y_versor,vector z_versor,sph_ws *Legendre_workspace,double *Legendre_associated_polynomials)
{
	interactionmap *ime=op->ime;

	calculateInteractionMapWithCutoffDistanceOrdered(cells,ime,pos,Box,range);

	// settiamo il massimo numero di neighbours
	int index;
	for (index=0;index<ime->num;index++)
	{
		//printf("%d\n",ime->howmany[index]);
		if (ime->howmany[index]>maxneighbours)
			ime->howmany[index]=maxneighbours;
	}

	vector v_12,xy_projection,z_projection;


	// scegliamo una terna di riferimento casuale in modo da poter partire anche dal cristallo

	//vector old_v_12;
	double inorm;
	double costeta,cosphi,sinphi,phi;
	int particle1,particle2,j,l,m;

	for (particle1=0;particle1<ime->num;particle1++)
	{

		for (j=0;j<ime->howmany[particle1];j++)
		{
			particle2=ime->with[particle1][j];

			// calcoliamo il versore

			v_12.x=pos[particle1].x-pos[particle2].x;
			v_12.y=pos[particle1].y-pos[particle2].y;
			v_12.z=pos[particle1].z-pos[particle2].z;

			v_12.x-=Box->x*rint(v_12.x/Box->x);
			v_12.y-=Box->y*rint(v_12.y/Box->y);
			v_12.z-=Box->z*rint(v_12.z/Box->z);

			inorm=1./sqrt(SQR(v_12.x)+SQR(v_12.y)+SQR(v_12.z));

			v_12.x*=inorm;
			v_12.y*=inorm;
			v_12.z*=inorm;

			costeta=SCALAR(&v_12,&z_versor);


			double z_projection_module=SCALAR(&v_12,&z_versor);

			z_projection.x=z_projection_module*z_versor.x;
			z_projection.y=z_projection_module*z_versor.y;
			z_projection.z=z_projection_module*z_versor.z;

			xy_projection.x=v_12.x-z_projection.x;
			xy_projection.y=v_12.y-z_projection.y;
			xy_projection.z=v_12.z-z_projection.z;

			inorm=1./sqrt(SQR(xy_projection.x)+SQR(xy_projection.y)+SQR(xy_projection.z));
			xy_projection.x*=inorm;
			xy_projection.y*=inorm;
			xy_projection.z*=inorm;

			cosphi=SCALAR(&xy_projection,&x_versor);
			sinphi=SCALAR(&xy_projection,&y_versor);

			if (sinphi>0)
				phi=acos(cosphi);
			else
				phi=2.*M_PI-acos(cosphi);


			complex double expphi=cexp(I*phi);
			complex double expumphi=cexp(-I*phi);

			int il;
			for (il=0;il<op->q_num;il++)
			{
				l=op->q_list[il];

				int order=ordinatorGetOrder(op->o,l);

				// precalcolo dei polinomi di legendre
				dPl0dm (Legendre_workspace,Legendre_associated_polynomials,costeta,l);


				op->q[order][particle1][l]+=Legendre_associated_polynomials[0]/(ime->howmany[particle1]);

				complex double exphi_last=expphi;
				complex double exumphi_last=expumphi;

				for (m=1;m<=l;m++)
				{
					complex double spherical_harmonic=Legendre_associated_polynomials[m]*exphi_last;

					op->q[order][particle1][m+l]+=spherical_harmonic/(ime->howmany[particle1]);

					spherical_harmonic=(m%2?-1:1)*Legendre_associated_polynomials[m]*exumphi_last;

					op->q[order][particle1][-m+l]+=spherical_harmonic/(ime->howmany[particle1]);

					exphi_last*=expphi;
					exumphi_last*=expumphi;
				}
			}

		}

	}

}

void calculate_Qlm(orderparam *op,int ncolloids)
{
// 	interactionmap *im=op->im;
	interactionmap *ime=op->ime;

	int il,particle1,particle2,m,j;
	for (il=0;il<op->Q_num;il++)
	{
		int l=op->Q_list[il];
		int order=ordinatorGetOrder(op->o,l);

		for (particle1=0;particle1<ncolloids;particle1++)
		{
			for (m=-l;m<=l;m++)
			{
				op->Q[order][particle1][m+l]+=op->q[order][particle1][m+l]/(1+ime->howmany[particle1]);
			}

			for (j=0;j<ime->howmany[particle1];j++)
			{
				particle2=ime->with[particle1][j];

				for (m=-l;m<=l;m++)
				{
					op->Q[order][particle1][m+l]+=op->q[order][particle2][m+l]/(1+ime->howmany[particle1]);
					//op->Q[order][particle2][m+l]+=op->q[order][particle1][m+l]/(1+ime->howmany[particle2]);
				}
			}
		}
	}
}


void calculate_QlmFlavio(orderparam *op,int ncolloids)
{

	interactionmap *ime=op->ime;

	int il,particle1,particle2,m,j;
	for (il=0;il<op->Q_num;il++)
	{
		int l=op->Q_list[il];
		int order=ordinatorGetOrder(op->o,l);

		for (particle1=0;particle1<ncolloids;particle1++)
		{
			for (m=-l;m<=l;m++)
			{
				op->Q[order][particle1][m+l]=0;
				op->Q[order][particle1][m+l]+=op->q[order][particle1][m+l]/(1+ime->howmany[particle1]);
			}

			for (j=0;j<ime->howmany[particle1];j++)
			{
				particle2=ime->with[particle1][j];

				for (m=-l;m<=l;m++)
				{


					op->Q[order][particle1][m+l]+=op->q[order][particle2][m+l]/(1+ime->howmany[particle1]);
					//op->Q[order][particle2][m+l]+=op->q[order][particle1][m+l]/(1+ime->howmany[particle2]);
				}
				for (m=-l;m<=l;m++)
				{
					printf("%d %d -> %lf %lf\n",j,particle2,creal(op->Q[order][particle1][m+l]),cimag(op->Q[order][particle1][m+l]));
				}
			}

			exit(1);
		}
	}
}

void calculate_qlm_different_im(orderparam *op,vector *pos,vector *Box,double range,listcell *cells,vector x_versor,vector y_versor,vector z_versor,interactionmap *ime)
{
// 	calculateExtendedInteractionMapWithCutoffDistance(cells,op->im,op->ime,pos,box,range);
// 	interactionmap *im=op->im;
// 	interactionmap *ime=op->ime;

	vector v_12,xy_projection,z_projection;
	//vector old_v_12;


	// scegliamo una terna di riferimento casuale in modo da poter partire anche dal cristallo


	double inorm;
	double costeta,cosphi,sinphi,phi;
	int particle1,particle2,j,l,m;

	for (particle1=0;particle1<ime->num;particle1++)
	{

		for (j=0;j<ime->howmany[particle1];j++)
		{
			particle2=ime->with[particle1][j];

			// calcoliamo il versore

			v_12.x=pos[particle1].x-pos[particle2].x;
			v_12.y=pos[particle1].y-pos[particle2].y;
			v_12.z=pos[particle1].z-pos[particle2].z;

			v_12.x-=Box->x*rint(v_12.x/Box->x);
			v_12.y-=Box->y*rint(v_12.y/Box->y);
			v_12.z-=Box->z*rint(v_12.z/Box->z);

			inorm=1./sqrt(SQR(v_12.x)+SQR(v_12.y)+SQR(v_12.z));

			v_12.x*=inorm;
			v_12.y*=inorm;
			v_12.z*=inorm;

			costeta=SCALAR(&v_12,&z_versor);


			double z_projection_module=SCALAR(&v_12,&z_versor);

			z_projection.x=z_projection_module*z_versor.x;
			z_projection.y=z_projection_module*z_versor.y;
			z_projection.z=z_projection_module*z_versor.z;

			xy_projection.x=v_12.x-z_projection.x;
			xy_projection.y=v_12.y-z_projection.y;
			xy_projection.z=v_12.z-z_projection.z;

			inorm=1./sqrt(SQR(xy_projection.x)+SQR(xy_projection.y)+SQR(xy_projection.z));
			xy_projection.x*=inorm;
			xy_projection.y*=inorm;
			xy_projection.z*=inorm;

			cosphi=SCALAR(&xy_projection,&x_versor);
			sinphi=SCALAR(&xy_projection,&y_versor);

			if (sinphi>0)
				phi=acos(cosphi);
			else
				phi=2.*M_PI-acos(cosphi);

			int il;
			for (il=0;il<op->q_num;il++)
			{
				l=op->q_list[il];

				int order=ordinatorGetOrder(op->o,l);

				// precalcolo dei polinomi di legendre
				dPl0dm (&Legendre_workspace,Legendre_associated_polynomials,costeta,l);

				for (m=-l;m<=l;m++)
				{

					complex double spherical_harmonic;

					/*
					if (m>=0)
					{
						spherical_harmonic=gsl_sf_legendre_sphPlm(l,m,costeta)*cexp(I*m*phi);
					}
					else
					{
						spherical_harmonic=pow(-1.,-m)*gsl_sf_legendre_sphPlm(l,-m,costeta)*cexp(I*m*phi);
					}
					*/
					if (m>=0)
					{
						spherical_harmonic=Legendre_associated_polynomials[m]*cexp(I*m*phi);
					}
					else
					{
						spherical_harmonic=((-m)%2?-1:1)*Legendre_associated_polynomials[-m]*cexp(I*m*phi);
					}


					op->q[order][particle1][m+l]+=spherical_harmonic/(ime->howmany[particle1]);
					//op->q[order][particle2][m+l]+=spherical_harmonic/(ime->howmany[particle2]);
					//op->q[order][particle2][m+l]+=pow(-1,l)*spherical_harmonic/(ime->howmany[particle2]);
				}
			}

		}

	}

}

void calculate_Qlm_different_im(orderparam *op,int ncolloids,interactionmap *ime)
{
// 	interactionmap *im=op->im;
// 	interactionmap *ime=op->ime;

	int il,particle1,particle2,m,j;
	for (il=0;il<op->Q_num;il++)
	{
		int l=op->Q_list[il];
		int order=ordinatorGetOrder(op->o,l);

		for (particle1=0;particle1<ncolloids;particle1++)
		{
			for (m=-l;m<=l;m++)
			{
				op->Q[order][particle1][m+l]+=op->q[order][particle1][m+l]/(1+ime->howmany[particle1]);
			}

			for (j=0;j<ime->howmany[particle1];j++)
			{
				particle2=ime->with[particle1][j];

				for (m=-l;m<=l;m++)
				{
					op->Q[order][particle1][m+l]+=op->q[order][particle2][m+l]/(1+ime->howmany[particle1]);
					//op->Q[order][particle2][m+l]+=op->q[order][particle1][m+l]/(1+ime->howmany[particle2]);
				}
			}
		}
	}
}

double* get_ql(orderparam *op,int l,int ncolloids)
{
	int particle1,m;

	int order=ordinatorGetOrder(op->o,l);

	for (particle1=0;particle1<ncolloids;particle1++)
	{
		double q=0.;

		for (m=-l;m<=l;m++)
		{
			q+=SQR(creal(op->q[order][particle1][m+l]))+SQR(cimag(op->q[order][particle1][m+l]));
		}

		q*=(4*M_PI/(2.*l+1));

		op->ql[order][particle1]=sqrt(q);

	}

	return op->ql[order];
}

double* get_Ql(orderparam *op,int l,int ncolloids)
{
	int particle1,m;

	int order=ordinatorGetOrder(op->o,l);

	for (particle1=0;particle1<ncolloids;particle1++)
	{
		double q=0.;

		for (m=-l;m<=l;m++)
		{
			q+=SQR(creal(op->Q[order][particle1][m+l]))+SQR(cimag(op->Q[order][particle1][m+l]));
		}

		q*=(4*M_PI/(2.*l+1));

		op->Ql[order][particle1]=sqrt(q);

	}

	return op->Ql[order];
}

complex double** get_Qlm(orderparam *op,int l)
{
	int order=ordinatorGetOrder(op->o,l);

	return op->Q[order];
}

double *get_wl(orderparam *op,int l,int ncolloids)
{
	double *norm=calloc(ncolloids,sizeof(double));


	int order=ordinatorGetOrder(op->o,l);

	int m1,m2,m3; // c'e' il vincolo m1+m2+m3=0

	for (m1=-l;m1<=l;m1++)
	{
		for (m2=-l;m2<=l;m2++)
		{
			m3=-m1-m2;

			if (abs(m3)<=l)
			{
				double wigner=gsl_sf_coupling_3j(2*l,2*l,2*l,2*m1,2*m2,2*m3);
				int i;
				for (i=0;i<ncolloids;i++)
				{
					double product=creal(op->q[order][i][m1+l]*op->q[order][i][m2+l]*op->q[order][i][m3+l]);

					op->wl[order][i]+=wigner*product;
				}
			}
		}
	}

	// normalizziamo
	int particle1,m;
	for (particle1=0;particle1<ncolloids;particle1++)
	{
		for (m=-l;m<=l;m++)
		{
			norm[particle1]+=SQR(creal(op->q[order][particle1][m+l]))+SQR(cimag(op->q[order][particle1][m+l]));
		}

		op->wl[order][particle1]/=pow(norm[particle1],3./2.);

		norm[particle1]=0.;
	}

	free(norm);

	return op->wl[order];
}

double *get_wl_notnormalized(orderparam *op,int l,int ncolloids)
{
	int order=ordinatorGetOrder(op->o,l);

	int m1,m2,m3; // c'e' il vincolo m1+m2+m3=0

	for (m1=-l;m1<=l;m1++)
	{
		for (m2=-l;m2<=l;m2++)
		{
			m3=-m1-m2;

			if (abs(m3)<=l)
			{
				double wigner=gsl_sf_coupling_3j(2*l,2*l,2*l,2*m1,2*m2,2*m3);
				int i;
				for (i=0;i<ncolloids;i++)
				{
					double product=creal(op->q[order][i][m1+l]*op->q[order][i][m2+l]*op->q[order][i][m3+l]);

					op->wl[order][i]+=wigner*product;
				}
			}
		}
	}

	return op->wl[order];
}


double *get_Wl(orderparam *op,int l,int ncolloids)
{
	double *norm=calloc(ncolloids,sizeof(double));

	int order=ordinatorGetOrder(op->o,l);

	int m1,m2,m3; // c'e' il vincolo m1+m2+m3=0

	for (m1=-l;m1<=l;m1++)
	{
		for (m2=-l;m2<=l;m2++)
		{
			m3=-m1-m2;

			if (abs(m3)<=l)
			{
				double wigner=gsl_sf_coupling_3j(2*l,2*l,2*l,2*m1,2*m2,2*m3);
				int i;
				for (i=0;i<ncolloids;i++)
				{
					double product=creal(op->Q[order][i][m1+l]*op->Q[order][i][m2+l]*op->Q[order][i][m3+l]);

					op->Wl[order][i]+=wigner*product;
				}
			}
		}
	}

	// normalizziamo
	int particle1,m;
	for (particle1=0;particle1<ncolloids;particle1++)
	{
		for (m=-l;m<=l;m++)
		{
			norm[particle1]+=SQR(creal(op->Q[order][particle1][m+l]))+SQR(cimag(op->Q[order][particle1][m+l]));
		}

		op->Wl[order][particle1]/=pow(norm[particle1],3./2.);

		norm[particle1]=0.;
	}

	free(norm);

	return op->Wl[order];
}

double *get_Wl_notnormalized(orderparam *op,int l,int ncolloids)
{
	int order=ordinatorGetOrder(op->o,l);

	int m1,m2,m3; // c'e' il vincolo m1+m2+m3=0

	for (m1=-l;m1<=l;m1++)
	{
		for (m2=-l;m2<=l;m2++)
		{
			m3=-m1-m2;

			if (abs(m3)<=l)
			{
				double wigner=gsl_sf_coupling_3j(2*l,2*l,2*l,2*m1,2*m2,2*m3);
				int i;
				for (i=0;i<ncolloids;i++)
				{
					double product=creal(op->Q[order][i][m1+l]*op->Q[order][i][m2+l]*op->Q[order][i][m3+l]);

					op->Wl[order][i]+=wigner*product;
				}
			}
		}
	}

	return op->Wl[order];
}

double get_ql_global(orderparam *op,int l,int ncolloids)
{
	interactionmap *ime=op->ime;
	int order=ordinatorGetOrder(op->o,l);

	double global_q=0.;
	int m;
	for (m=-l;m<=l;m++)
	{
		complex double Qlm=0.;
		int totbonds=0;

		int i;
		for (i=0;i<ncolloids;i++)
		{
			Qlm+=op->q[order][i][m+l]*(ime->howmany[i]);
			totbonds+=ime->howmany[i];
		}

		Qlm/=(double)totbonds;

		global_q+=SQR(creal(Qlm))+SQR(cimag(Qlm));
	}

	//global_q*=(4*M_PI/(2*l+1));


	return sqrt(global_q);
}

double get_Ql_global(orderparam *op,int l,int ncolloids)
{
	interactionmap *ime=op->ime;
	int order=ordinatorGetOrder(op->o,l);

	double global_Q=0.;
	int m;
	for (m=-l;m<=l;m++)
	{
		complex double Qlm=0.;
		int totbonds=0;

		int i;
		for (i=0;i<ncolloids;i++)
		{
			Qlm+=op->Q[order][i][m+l]*(ime->howmany[i]);
			totbonds+=ime->howmany[i];
		}

		Qlm/=(double)totbonds;

		global_Q+=SQR(creal(Qlm))+SQR(cimag(Qlm));
	}

	//global_Q*=(4*M_PI/(2*l+1));


	return sqrt(global_Q);
}

void getCoherentParticles(int l,complex double ***q_local,orderparam *op,int ncolloids,double coherence_threshold,interactionmap *coherent_map,int *num_coherent,double *buffer)
{
	// shortcuts
	double *q_local_norm=buffer;
	interactionmap *ime=op->ime;

	int order=ordinatorGetOrder(op->o,l);

	int i,m,particle1,particle2,j;

	resetInteractionMap(coherent_map);

	// normalizziamo la norma di q_local
	for (i=0;i<ncolloids;i++)
	{
		q_local_norm[i]=0.;

		for (m=-l;m<=l;m++)
		{
			q_local_norm[i]+=SQR(creal(q_local[order][i][m+l]))+SQR(cimag(q_local[order][i][m+l]));
		}
	}


	// troviamo le coppie di particelle coerenti

	for (i=0;i<ncolloids;i++)
	{
		num_coherent[i]=0;
	}

	coherent_map->num=ncolloids;


	for (particle1=0;particle1<ime->num;particle1++)
	{
		int howmany_particle1=0;

		for (j=0;j<ime->howmany[particle1];j++)
		{
			particle2=ime->with[particle1][j];

			// prodotto scalare

			complex double sp=0.+0.*I;

			for (m=-l;m<=l;m++)
			{
				sp+=q_local[order][particle1][m+l]*conj(q_local[order][particle2][m+l]);
			}

			double coherence=creal(sp)/sqrt(q_local_norm[particle1]*q_local_norm[particle2]);


			if (coherence>coherence_threshold)
			{
				num_coherent[particle1]++;
				//num_coherent[particle2]++;

				// si costruisce una lista di particelle interagenti
				coherent_map->who[particle1]=particle1;
				coherent_map->with[particle1][howmany_particle1++]=particle2;
			}
		}

		coherent_map->howmany[particle1]=howmany_particle1;
	}
}

void getCoherentParticlesAndTotalCoherence(int l,complex double ***q_local,orderparam *op,int ncolloids,double coherence_threshold,interactionmap *coherent_map,int *num_coherent,double *total_coherence,double *buffer)
{
	// shortcuts
	double *q_local_norm=buffer;
	interactionmap *ime=op->ime;

	int order=ordinatorGetOrder(op->o,l);

	int i,m,particle1,particle2,j;

	resetInteractionMap(coherent_map);

	// normalizziamo la norma di q_local
	for (i=0;i<ncolloids;i++)
	{
		q_local_norm[i]=0.;

		for (m=-l;m<=l;m++)
		{
			q_local_norm[i]+=SQR(creal(q_local[order][i][m+l]))+SQR(cimag(q_local[order][i][m+l]));
		}
	}


	// troviamo le coppie di particelle coerenti

	for (i=0;i<ncolloids;i++)
	{
		num_coherent[i]=0;
	}

	coherent_map->num=ncolloids;


	for (particle1=0;particle1<ime->num;particle1++)
	{
		int howmany_particle1=0;

		total_coherence[particle1]=0.;

		for (j=0;j<ime->howmany[particle1];j++)
		{
			particle2=ime->with[particle1][j];

			// prodotto scalare

			complex double sp=0.+0.*I;

			for (m=-l;m<=l;m++)
			{
				sp+=q_local[order][particle1][m+l]*conj(q_local[order][particle2][m+l]);
			}

			double coherence=creal(sp)/sqrt(q_local_norm[particle1]*q_local_norm[particle2]);

			total_coherence[particle1]+=coherence;

			if (coherence>coherence_threshold)
			{
				num_coherent[particle1]++;
				//num_coherent[particle2]++;

				// si costruisce una lista di particelle interagenti
				coherent_map->who[particle1]=particle1;
				coherent_map->with[particle1][howmany_particle1++]=particle2;
			}
		}

		coherent_map->howmany[particle1]=howmany_particle1;
	}
}


int request_OpL(char op_code,int l,int ncolloids,orderparam *op)
{
	if (op_code=='q')
		return request_ql(l,ncolloids,op);
	else if (op_code=='Q')
		return request_Ql(l,ncolloids,op);
	else if (op_code=='w')
		return request_wl(l,ncolloids,op);
	else if (op_code=='W')
		return request_Wl(l,ncolloids,op);
	else
	{
		printf("Error: wrong order parameter code\n");
		exit(1);
		return 0;
	}
}

double* get_OpL(char op_code,int l,int ncolloids,orderparam *op)
{
	if (op_code=='q')
		return get_ql(op,l,ncolloids);
	else if (op_code=='Q')
		return get_Ql(op,l,ncolloids);
	else if (op_code=='w')
		return get_wl_notnormalized(op,l,ncolloids);
	else if (op_code=='W')
		return get_Wl(op,l,ncolloids);
	else
	{
		printf("Error: wrong order parameter code\n");
		exit(1);
		return 0;
	}
}

double get_OpL_global(char op_code,int l,int ncolloids,orderparam *op)
{
	if (op_code=='q')
		return get_ql_global(op,l,ncolloids);
	else if (op_code=='Q')
		return get_Ql_global(op,l,ncolloids);
	else
	{
		printf("Error: wrong order parameter code\n");
		exit(1);
		return 0;
	}
}


/* PUBLIC INTERFACE */
void crystalsConstructor(input_file *input,int ncolloids,vector box,Output *output_files)
{
    Method=0;

    int found_method=getInputInt("Order_parameter_method",&Method,0);

    if ((found_method==KEY_NOT_FOUND) || (Method==0))
    {
        output_log_msg(,output_files,"Order Parameter calculation disabled\n");
        return;
    }

	getInputDouble("Order_parameter_neighbours_range",&Range,1);
	getInputInt("Order_parameter_max_neighbours",&MaxNeighbours,1);
	
	getInputInt("Order_parameter_l",&OP_l,1);
	getInputString("Order_parameter_code",&OP_code,1);


	int found_threshold=getInputDouble("Order_parameter_coherence_threshold",&Coherence_threshold,0);
	int found_solid_threshold=getInputInt("Order_parameter_solid_threshold",&Solid_threshold,0);
	int found_staggered_threshold=getInputDouble("Order_parameter_staggered_threshold",&Staggered_threshold,0);
	int found_eclipsed_lowthreshold=getInputDouble("Order_parameter_eclipsed_lowthreshold",&Eclipsed_lowthreshold,0);
	int found_eclipsed_highthreshold=getInputDouble("Order_parameter_eclipsed_highthreshold",&Eclipsed_highthreshold,0);
	

	CrystalFractionCell=getList(box,Range,ncolloids);

    OP=allocateOP(ncolloids,3);

	request_OpL(OP_code,OP_l,ncolloids,OP);

	if (OP_code=='q')
		Punt_qlm=OP->q;
	else if (OP_code=='Q')
		Punt_qlm=OP->Q;
	else
	{
		logPrint("Error: Order_parameter_code can be either 'q' or 'Q'\n");
		exit(1);
	}

	sph_ws_init(&Legendre_workspace,OP_l);

	Legendre_associated_polynomials=calloc(OP_l+1,sizeof(double));


	// costruiamo il sistema di riferimento
    random_vector_on_sphere(X_versor);
    get_perpendicular_versor(X_versor,Y_versor);
	cross(X_versor,Y_versor,Z_versor);

	Buffer=calloc(ncolloids,sizeof(double));


    output_log_msg(output_files,"\nORDER PARAMETER SECTION\n");
	output_log_msg(output_files,"Order parameter %c%d\n",OP_code,OP_l);
	output_log_msg(output_files,"Neighour's range %lf\n",Range);
	output_log_msg(output_files,"Max neighbours %d\n",MaxNeighbours);

	ImSolid=calloc(ncolloids,sizeof(int));
	Cluster_solid=getClusters(ncolloids);
	List_solid=bilistaGet(ncolloids);

	
    // FRENKEL-like METHOD

    if ((found_threshold==KEY_NOT_FOUND) || (found_solid_threshold==KEY_NOT_FOUND))
    {
        output_log_msg(,output_files,"Error: define 'Order_parameter_coherence_threshold' and 'Order_parameter_solid_threshold'\n");
        exit(1);
    }
    else
    {
        output_log_msg("Metodo: Frenkel\n");
        output_log_msg("Coherence threshold %lf\n",Coherence_threshold);
        output_log_msg("Number links per solid particle %d\n\n",Solid_threshold);
    }

    Coherent_map=createInteractionMap(ncolloids,MaxNeighbours);
    Num_coherent=calloc(ncolloids,sizeof(int));


	searchFree(s);
}


void freeCrystals()
{
	if (Method==0)
		return;

	freeOP(OP,3);

	free(ImSolid);
	freeClusters(Cluster_solid);
	bilistaFree(List_solid);

	freeInteractionMap(Coherent_map);
	free(Num_coherent);
	free(Buffer);

	freeList(CrystalFractionCell);

	sph_ws_free(&Legendre_workspace);
	free(Legendre_associated_polynomials);
}


double getOrderParameter(System *syst,int *num_solid)
{
    // vector *pos,int ncolloids,vector *box


	int particle1,particle2,j;


	if (Method==-1)
	{
		*num_solid=0.;
		return -1;
	}

	resetOP(OP,ncolloids,2);


	fullUpdateList(CrystalFractionCell,pos,ncolloids,*box,Range);

	calculate_qlm_maxneighbours_fast(OP,pos,box,Range,MaxNeighbours,CrystalFractionCell,X_versor,Y_versor,Z_versor,&Legendre_workspace,Legendre_associated_polynomials);
	
    calculate_Qlm(OP,ncolloids);


	
    getCoherentParticles(OP_l,Punt_qlm,OP,ncolloids,Coherence_threshold,Coherent_map,Num_coherent,Buffer);

    for (particle1=0;particle1<ncolloids;particle1++)
    {

        if (Num_coherent[particle1]>=Solid_threshold)
        {
            ImSolid[particle1]=1;
        }
        else
        {
            ImSolid[particle1]=0;
        }

    }

	
	

	resetClusters(Cluster_solid,ncolloids);
    bilistaReset(List_solid,ncolloids);

	*num_solid=0;
	int max_size=0;
	for (particle1=0;particle1<ncolloids;particle1++)
	{

		if (ImSolid[particle1]==1)
		{
			bilistaInsert(List_solid,particle1);
			(*num_solid)++;

			addNode(particle1,Cluster_solid);

			if (Method==0)
			{
				// USUAL
				interactionmap *ime=OP->ime;
				for (j=0;j<ime->howmany[particle1];j++)
				{
					particle2=ime->with[particle1][j];

				// SAIKA
// 				for (j=0;j<Coherent_map->howmany[particle1];j++)
// 				{
// 					particle2=Coherent_map->with[particle1][j];

					if (ImSolid[particle2]==1)
					{

						addNode(particle2,Cluster_solid);

						int size=addBond(particle1,particle2,Cluster_solid);

						if (size>max_size)
							max_size=size;

					}
				}
			}
			else if ((Method==1) || (Method==2))
			{
				for (j=0;j<Staggered_map->howmany[particle1];j++)
				{
					particle2=Staggered_map->with[particle1][j];

					if (ImSolid[particle2]==1)
					{

						addNode(particle2,Cluster_solid);

						int size=addBond(particle1,particle2,Cluster_solid);

						if (size>max_size)
							max_size=size;

					}
				}
				for (j=0;j<Eclipsed_map->howmany[particle1];j++)
				{
					particle2=Eclipsed_map->with[particle1][j];

					if (ImSolid[particle2]==1)
					{

						addNode(particle2,Cluster_solid);

						int size=addBond(particle1,particle2,Cluster_solid);

						if (size>max_size)
							max_size=size;

					}
				}
			}
			else if (Method==3)
			{
				interactionmap *ime=OP->ime;

				for (j=0;j<ime->howmany[particle1];j++)
				{
					particle2=ime->with[particle1][j];

					if (ImSolid[particle2]==1)
					{

						addNode(particle2,Cluster_solid);

						int size=addBond(particle1,particle2,Cluster_solid);

						if (size>max_size)
							max_size=size;

					}
				}
			}
		}
	}

	if ( (max_size==0) && (*num_solid>0) )
		max_size=1;

	return (double)max_size;
}

