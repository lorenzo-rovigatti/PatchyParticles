#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


static const gsl_rng_type * Gsl_T;
static gsl_rng * Gsl_random;


typedef double vector[3];
typedef double matrix[3][3];


#define SQR(x) ((x) * (x))
#define SCALAR(x, y) ((x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2])

#define MATRIX_VECTOR_MULTIPLICATION(m, v, result) {\
	(result)[0] = (m)[0][0]*(v)[0] + (m)[0][1]*(v)[1] + (m)[0][2]*(v)[2];\
	(result)[1] = (m)[1][0]*(v)[0] + (m)[1][1]*(v)[1] + (m)[1][2]*(v)[2];\
	(result)[2] = (m)[2][0]*(v)[0] + (m)[2][1]*(v)[1] + (m)[2][2]*(v)[2];\
}


unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

void normalize(vector v) {
	double norm = sqrt(SCALAR(v, v));
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}


void set_vector(vector v, double x, double y, double z) {
	v[0] = x;
	v[1] = y;
	v[2] = z;
}


double determinant(double (*m)[3]) {
	double det = 0.;

	det += (m[0][0])*(m[1][1])*(m[2][2]);
	det += (m[0][1])*(m[1][2])*(m[2][0]);
	det += (m[0][2])*(m[1][0])*(m[2][1]);
	det -= (m[0][2])*(m[1][1])*(m[2][0]);
	det -= (m[0][0])*(m[1][2])*(m[2][1]);
	det -= (m[0][1])*(m[1][0])*(m[2][2]);

	return det;
}

void get_rotation_matrix(vector axis, double t, matrix rotation_matrix) {
	double st = sin(t);
	double ct = cos(t);
	double uct = 1. - ct;

	rotation_matrix[0][0] = SQR(axis[0]) + (1. - SQR(axis[0]))*ct;
	rotation_matrix[0][1] = axis[0]*axis[1]*uct - axis[2]*st;
	rotation_matrix[0][2] = axis[0]*axis[2]*uct + axis[1]*st;

	rotation_matrix[1][0] = axis[0]*axis[1]*uct + axis[2]*st;
	rotation_matrix[1][1] = SQR(axis[1]) + (1. - SQR(axis[1]))*ct;
	rotation_matrix[1][2] = axis[1]*axis[2]*uct - axis[0]*st;

	rotation_matrix[2][0] = axis[0]*axis[2]*uct - axis[1]*st;
	rotation_matrix[2][1] = axis[1]*axis[2]*uct + axis[0]*st;
	rotation_matrix[2][2] = SQR(axis[2]) + (1. - SQR(axis[2]))*ct;
}

void get_rotation_matrix_2(vector axis, double t, matrix rotation_matrix) {
	double w = cos(0.5*t);
	double x = axis[0]*sin(0.5*t);
  double y = axis[1]*sin(0.5*t);
  double z = axis[2]*sin(0.5*t);

	rotation_matrix[0][0] = 1-2*(y*y+z*z);
	rotation_matrix[0][1] = 2*x*y+2*w*z;
	rotation_matrix[0][2] = 2*x*z-2*w*y;

	rotation_matrix[1][0] = 2*x*y-2*w*z;
	rotation_matrix[1][1] = 1-2*(x*x+z*z);
	rotation_matrix[1][2] = 2*y*z+2*w*x;

	rotation_matrix[2][0] = 2*x*z+2*w*y;
	rotation_matrix[2][1] = 2*y*z-2*w*x;
	rotation_matrix[2][2] = 1-2*(y*y+x*x);
}


void random_vector_on_sphere(vector res) {
	double ransq;
	double ran1, ran2;
	double ranh;

	do {
		//ran1 = 1. - 2.*drand48();
		//ran2 = 1. - 2.*drand48();
    ran1 = 1. - 2.*gsl_rng_uniform(Gsl_random);
		ran2 = 1. - 2.*gsl_rng_uniform(Gsl_random);
		ransq = SQR(ran1) + SQR(ran2);
	} while(ransq >= 1.);

	ranh = 2.*sqrt(1. - ransq);

	res[0] = ran1*ranh;
	res[1] = ran2*ranh;
	res[2] = 1. - 2.*ransq;
}


void gram_schmidt(vector v1, vector v2, vector v3) {
	double v1_norm2 = SCALAR(v1, v1);
	double v2_v1 = SCALAR(v2, v1);
	double buffer1 = v2_v1 / v1_norm2;

	v2[0] -= buffer1 * v1[0];
	v2[1] -= buffer1 * v1[1];
	v2[2] -= buffer1 * v1[2];

	double v3_v1 = SCALAR(v3, v1);
	double v3_v2 = SCALAR(v3, v2);
	double v2_norm2 = SCALAR(v2, v2);
	buffer1 = v3_v1 / v1_norm2;
	double buffer2 = v3_v2 / v2_norm2;

	v3[0] -= buffer1*v1[0] + buffer2*v2[0];
	v3[1] -= buffer1*v1[1] + buffer2*v2[1];
	v3[2] -= buffer1*v1[2] + buffer2*v2[2];

	normalize(v1);
	normalize(v2);
	normalize(v3);
}


void random_orientation(matrix orient) {
	vector axis;
	random_vector_on_sphere(axis);
	//double t = drand48() * 2 * M_PI;
  double t = gsl_rng_uniform(Gsl_random) * 2 * M_PI;
	get_rotation_matrix(axis, t, orient);
}


void random_orientation_2(matrix orient) {

	vector x,y,z;

  random_vector_on_sphere(x);
  random_vector_on_sphere(y);
  random_vector_on_sphere(z);

  // orthonormalize
  gram_schmidt(x,y,z);



  orient[0][0] = z[0];
  orient[0][1] = z[1];
  orient[0][2] = z[2];

  orient[1][0] = x[0];
  orient[1][1] = x[1];
  orient[1][2] = x[2];

  orient[2][0] = y[0];
  orient[2][1] = y[1];
  orient[2][2] = y[2];

  // rotations have det(R) == 1
  if(determinant(orient) < 0) {
          orient[0][0] = x[0];
          orient[0][1] = x[1];
          orient[0][2] = x[2];

          orient[1][0] = z[0];
          orient[1][1] = z[1];
          orient[1][2] = z[2];
  }

}


void random_orientation_3(matrix rotation_matrix) {
  double w = gsl_ran_gaussian(Gsl_random,1.);
	double x = gsl_ran_gaussian(Gsl_random,1.);
  double y = gsl_ran_gaussian(Gsl_random,1.);
  double z = gsl_ran_gaussian(Gsl_random,1.);

  double norm=sqrt(SQR(w)+SQR(x)+SQR(y)+SQR(z));

  w/=norm;
  x/=norm;
  y/=norm;
  z/=norm;

	rotation_matrix[0][0] = 1-2*(y*y+z*z);
	rotation_matrix[0][1] = 2*x*y+2*w*z;
	rotation_matrix[0][2] = 2*x*z-2*w*y;

	rotation_matrix[1][0] = 2*x*y-2*w*z;
	rotation_matrix[1][1] = 1-2*(x*x+z*z);
	rotation_matrix[1][2] = 2*y*z+2*w*x;

	rotation_matrix[2][0] = 2*x*z+2*w*y;
	rotation_matrix[2][1] = 2*y*z-2*w*x;
	rotation_matrix[2][2] = 1-2*(y*y+x*x);

}


int main(int argc,char *argv[])
{

  gsl_rng_env_setup();
  Gsl_T = gsl_rng_rand48;
  Gsl_random = gsl_rng_alloc(Gsl_T);

  int i;

  //unsigned long seed = mix(clock(), time(NULL), getpid());
  //srand48(seed);

  vector *base_patches = malloc(sizeof(vector) * 4);
  vector *real_patches = malloc(sizeof(vector) * 4);

	double half_isqrt3 = 0.5 / sqrt(3);
	set_vector(base_patches[0], -half_isqrt3, -half_isqrt3,  half_isqrt3);
	set_vector(base_patches[1], half_isqrt3, -half_isqrt3, -half_isqrt3);
	set_vector(base_patches[2], half_isqrt3,  half_isqrt3,  half_isqrt3);
	set_vector(base_patches[3], -half_isqrt3,  half_isqrt3, -half_isqrt3);
  for(i = 0; i < 4; i++) normalize(base_patches[i]);

  set_vector(base_patches[1], 0., 0.,  1.);

  matrix orientazione;

  int j;
  int sum=0;
  int sum2=0;
  for (i=0;i<10000000;i++)
  {
    random_orientation_2(orientazione);

    for(j=0;j<4;j++) MATRIX_VECTOR_MULTIPLICATION(orientazione,base_patches[j],real_patches[j]);

    //if (real_patches[0][2]>0.999)
    //  sum2+=1;

    //if (real_patches[1][0]*base_patches[0][0]+real_patches[1][1]*base_patches[0][1]+real_patches[1][2]*base_patches[0][2]>0.999)
    //  sum+=1;

    printf("%lf\n",real_patches[0][2]);

  }

  //printf("%d %d\n",sum,sum2);

  return 0;
}
