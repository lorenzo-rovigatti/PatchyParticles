#ifndef ORDER_PARAMETERS_H
#define ORDER_PARAMETERS_H

#define JSCALAR(v,u) (((v)->x)*((u)->x)+((v)->y)*((u)->y)+((v)->z)*((u)->z))

typedef struct _orderparam {
	ordinator *o;

	interactionmap *im;
	interactionmap *ime;

	complex double ***q;
	int *q_list;
	int q_num;

	complex double ***Q;
	int *Q_list;
	int Q_num;

	double **ql;

	double **Ql;

	double **wl;

	double **Wl;
} orderparam;


// per il calcolo veloce dei polinomi di Legendre
typedef struct {
	int lmax;
	double **PL0, **nrm;
	double *tmp, *xp;
} sph_ws;

// inline double intpow (const double x, const int i);
// inline double cintpow (const complex double x, const int i);


void crystalsConstructor(input_file *input,Output *output_files,System *syst);
void freeCrystals();
double getOrderParameter(System *syst,int *num_solid);



#endif