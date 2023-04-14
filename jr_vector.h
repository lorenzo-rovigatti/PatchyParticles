#ifndef VECTOR_H
#define VECTOR_H

typedef struct _vector {
	double x;
	double y;
	double z;
} vector;

vector vectorAdd(const vector *v1,const vector *v2);
vector vectorSub(const vector *v1,const vector *v2);
vector* vectorScale (double scalar,vector *v);
double vectorScalarProduct (const vector *v1,const vector *v2);
double vectorNorm (const vector *v);
double vectorSquareNorm (const vector *v);
vector* vectorVersor (vector *v);
vector* vectorOpposite(vector *v);
vector vectorVectorProduct (const vector *v1,const vector *v2);
void gramSchmidt(vector *v1,vector *v2,vector *v3);
double determinant(double (*m)[3]);
vector matrix_vector_multiplication(double (*m)[3],vector *v);
void randomVector(vector *rv,gsl_rng *random);
vector rotateVector(vector *v,vector *axis,double teta);
#endif

