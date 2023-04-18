#ifndef VECTOR_H
#define VECTOR_H

typedef struct _jvector {
	double x;
	double y;
	double z;
} jvector;

jvector vectorAdd(const jvector *v1,const jvector *v2);
jvector vectorSub(const jvector *v1,const jvector *v2);
jvector* vectorScale (double scalar,jvector *v);
double vectorScalarProduct (const jvector *v1,const jvector *v2);
double vectorNorm (const jvector *v);
double vectorSquareNorm (const jvector *v);
jvector* vectorVersor (jvector *v);
jvector* vectorOpposite(jvector *v);
jvector vectorVectorProduct (const jvector *v1,const jvector *v2);
void gramSchmidt(jvector *v1,jvector *v2,jvector *v3);
double jr_determinant(double (*m)[3]);
jvector matrix_vector_multiplication(double (*m)[3],jvector *v);
jvector rotateVector(jvector *v,jvector *axis,double teta);

#endif

