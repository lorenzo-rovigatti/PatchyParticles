/*
 * utils.h
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#ifndef UTILS_H_
#define UTILS_H_

#define COS120 -0.5
#define SIN120 0.8660254037844386

#define SQR(x) ((x) * (x))
#define SCALAR(x, y) ((x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2])

#define MATRIX_VECTOR_MULTIPLICATION(m, v, result) {\
	(result)[0] = (m)[0][0]*(v)[0] + (m)[0][1]*(v)[1] + (m)[0][2]*(v)[2];\
	(result)[1] = (m)[1][0]*(v)[0] + (m)[1][1]*(v)[1] + (m)[1][2]*(v)[2];\
	(result)[2] = (m)[2][0]*(v)[0] + (m)[2][1]*(v)[1] + (m)[2][2]*(v)[2];\
}

#include "defs.h"

void set_patches(System *syst, PatchyParticle *p);

void set_vector(vector v, double x, double y, double z);
void set_base_orientation(matrix orient);
void cross(vector v1, vector v2, vector res);
void normalize(vector v);
void matrix_matrix_multiplication(matrix m, matrix n, matrix res);
//void matrix_vector_multiplication(matrix m, vector v, vector result);
void random_vector_on_sphere(vector res);
void random_orientation(System *syst, matrix orient);
void get_rotation_matrix(vector axis, double t, matrix rotation_matrix);
void place_inside_vbonding(System *syst, PatchyParticle *rec, vector r, matrix orient, int rec_patch);
void utils_rotate_matrix(matrix orient_old, matrix orient_new, vector axis, double t);
void get_rotated_vector(vector v, vector axis, double t, vector res);
void rotate_vector(vector v, vector axis, double t);
void set_orientation_around_vector(vector v, matrix orient, double t);
void gram_schmidt(vector v1, vector v2, vector v3);
void utils_reset_acceptance_counters(System *syst);

#endif /* UTILS_H_ */
