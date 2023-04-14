#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "vector.h"

#define SQR(x) ((x)*(x))

vector vectorAdd(const vector *v1,const vector *v2)
{
	vector sum;
	
	sum.x=v1->x+v2->x;
	sum.y=v1->y+v2->y;
	sum.z=v1->z+v2->z;
	
	return sum;
}

vector vectorSub(const vector *v1,const vector *v2)
{
	vector sub;
	
	sub.x=v1->x-v2->x;
	sub.y=v1->y-v2->y;
	sub.z=v1->z-v2->z;
	
	return sub;
}

vector* vectorScale (double scalar,vector *v)
{
	v->x=scalar*v->x;
	v->y=scalar*v->y;
	v->z=scalar*v->z;
	
	return v;
}

double vectorScalarProduct (const vector *v1,const vector *v2)
{
	return v1->x*v2->x+v1->y*v2->y+v1->z*v2->z;
}

double vectorNorm (const vector *v)
{
	return sqrt(SQR(v->x)+SQR(v->y)+SQR(v->z));
}

double vectorSquareNorm (const vector *v)
{	
	return SQR(v->x)+SQR(v->y)+SQR(v->z);
}

vector* vectorVersor (vector *v)
{
	
	double i_n=1./vectorNorm(v);
	
	v->x=i_n*v->x;
	v->y=i_n*v->y;
	v->z=i_n*v->z;
	
	return v;
}

vector* vectorOpposite(vector *v)
{
	v->x=-v->x;
	v->y=-v->y;
	v->z=-v->z;
	return v;
}

vector vectorVectorProduct (const vector *v1,const vector *v2)
{
	vector result;
	
	result.x=(v1->y*v2->z)-(v1->z*v2->y);
	result.y=(v1->z*v2->x)-(v1->x*v2->z);
	result.z=(v1->x*v2->y)-(v1->y*v2->x);
	
	return result;
	
}

void gramSchmidt(vector *v1,vector *v2,vector *v3)
{
	double v1_norm2=SQR(v1->x)+SQR(v1->y)+SQR(v1->z);
	double v2_v1=v2->x*v1->x+v2->y*v1->y+v2->z*v1->z;
	double buffer=v2_v1/v1_norm2;
	
	
	// u2=v2-((v2*u1)/(u1.norma**2))*u1
	v2->x=v2->x-(buffer)*v1->x;
	v2->y=v2->y-(buffer)*v1->y;
	v2->z=v2->z-(buffer)*v1->z;
	
	double v3_v1=v3->x*v1->x+v3->y*v1->y+v3->z*v1->z;
	double v3_v2=v3->x*v2->x+v3->y*v2->y+v3->z*v2->z;
	double v2_norm2=SQR(v2->x)+SQR(v2->y)+SQR(v2->z);
	double buffer1=v3_v1/v1_norm2;
	double buffer2=v3_v2/v2_norm2;
	
	// u3=v3-((v3*u1)/(u1.norma**2))*u1-((v3*u2)/(u2.norma**2))*u2
	v3->x=v3->x-buffer1*v1->x-buffer2*v2->x;
	v3->y=v3->y-buffer1*v1->y-buffer2*v2->y;
	v3->z=v3->z-buffer1*v1->z-buffer2*v2->z;
	
}

double determinant(double (*m)[3])
{
	double det=0.;
	
	det=det+(m[0][0])*(m[1][1])*(m[2][2]);
	det=det+(m[0][1])*(m[1][2])*(m[2][0]);
	det=det+(m[0][2])*(m[1][0])*(m[2][1]);
	det=det-(m[0][2])*(m[1][1])*(m[2][0]);
	det=det-(m[0][0])*(m[1][2])*(m[2][1]);
	det=det-(m[0][1])*(m[1][0])*(m[2][2]);
	
	return det;
}

vector matrix_vector_multiplication(double (*m)[3],vector *v)
{
	vector result;
	
	result.x=m[0][0]*v->x+m[0][1]*v->y+m[0][2]*v->z;
	result.y=m[1][0]*v->x+m[1][1]*v->y+m[1][2]*v->z;
	result.z=m[2][0]*v->x+m[2][1]*v->y+m[2][2]*v->z;
	
	return result;
}

// Frenkel-Smit algorithm
void randomVector(vector *rv,gsl_rng *random)
{
	double ransq=2.;
	double ran1,ran2;
	double ranh;
	
	while (ransq>=1.)
	{
		ran1=1.-2.*gsl_rng_uniform(random);
		ran2=1.-2.*gsl_rng_uniform(random);
		ransq=SQR(ran1)+SQR(ran2);
	}
	
	ranh=2.*sqrt(1.-ransq);
	
	rv->x=ran1*ranh;
	rv->y=ran2*ranh;
	rv->z=1.-2.*ransq;
	
}

vector rotateVector(vector *v,vector *axis,double teta)
{
	double uct=1.-cos(teta);
	double st=sin(teta);
	double ct=cos(teta);
	
	vector vettorei=*axis;
	vectorVersor(&vettorei);
	
	double rotation_matrix[3][3];
	
	rotation_matrix[0][0]=SQR(vettorei.x)+(1.-SQR(vettorei.x))*ct;
	rotation_matrix[0][1]=vettorei.x*vettorei.y*uct-vettorei.z*st;
	rotation_matrix[0][2]=vettorei.x*vettorei.z*uct+vettorei.y*st;
	
	rotation_matrix[1][0]=vettorei.x*vettorei.y*uct+vettorei.z*st;
	rotation_matrix[1][1]=SQR(vettorei.y)+(1.-SQR(vettorei.y))*ct;
	rotation_matrix[1][2]=vettorei.y*vettorei.z*uct-vettorei.x*st;
	
	rotation_matrix[2][0]=vettorei.x*vettorei.z*uct-vettorei.y*st;
	rotation_matrix[2][1]=vettorei.y*vettorei.z*uct+vettorei.x*st;
	rotation_matrix[2][2]=SQR(vettorei.z)+(1.-SQR(vettorei.z))*ct;
	
	return matrix_vector_multiplication(rotation_matrix,v);
}
