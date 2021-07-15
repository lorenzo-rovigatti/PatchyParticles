/*
 * utils.c
 *
 *  Created on: 01/nov/2011
 *      Author: lorenzo
 */

#include "utils.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void set_vector(vector v, double x, double y, double z) {
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

void set_patches(System *syst, PatchyParticle *p) {
	int i;
	for(i = 0; i < syst->n_patches; i++) MATRIX_VECTOR_MULTIPLICATION(p->orientation, syst->base_patches[i], p->patches[i]);
}

void cross(vector v1, vector v2, vector res) {
	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void normalize(vector v) {
	double norm = sqrt(SCALAR(v, v));
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}

void get_perpendicular_versor(vector v, vector res) {
	random_vector_on_sphere(res);

	double v_norm2 = SCALAR(v, v);
	double res_v = SCALAR(res, v);
	double buffer = res_v / v_norm2;

	res[0] -= buffer * v[0];
	res[1] -= buffer * v[1];
	res[2] -= buffer * v[2];

	normalize(res);
}

void place_inside_vbonding(System *syst, PatchyParticle *rec, vector r, matrix orient, int target_patch) {
	utils_rotate_matrix(syst->base_orient, orient, syst->base_patches[0], 2*drand48()*M_PI);

	vector target_patch_dir;
	memcpy(target_patch_dir, rec->patches[target_patch], 3*sizeof(double));

	vector norm_vect;
	get_perpendicular_versor(target_patch_dir, norm_vect);

	// choose a new direction
	double theta = acos(syst->kf_cosmax + (1. - syst->kf_cosmax) * drand48());
	rotate_vector(target_patch_dir, norm_vect, theta);
	normalize(target_patch_dir);

	// another random angle
	double theta2 = acos(syst->kf_cosmax + (1. - syst->kf_cosmax) * drand48());
	set_orientation_around_vector(target_patch_dir, orient, theta2);

	// choose the distance
	double dist = pow(1. + drand48()*(syst->kf_delta*syst->kf_delta*syst->kf_delta + 3.*SQR(syst->kf_delta) + 3.*syst->kf_delta), 1. / 3.);

	r[0] = target_patch_dir[0]*dist + rec->r[0];
	r[1] = target_patch_dir[1]*dist + rec->r[1];
	r[2] = target_patch_dir[2]*dist + rec->r[2];
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

void random_vector_on_sphere(vector res) {
	double ransq;
	double ran1, ran2;
	double ranh;

	do {
		ran1 = 1. - 2.*drand48();
		ran2 = 1. - 2.*drand48();
		ransq = SQR(ran1) + SQR(ran2);
	} while(ransq >= 1.);

	ranh = 2.*sqrt(1. - ransq);

	res[0] = ran1*ranh;
	res[1] = ran2*ranh;
	res[2] = 1. - 2.*ransq;
}

void random_orientation(System *syst, matrix orient) {
	vector axis;
	random_vector_on_sphere(axis);
	double t = drand48() * 2 * M_PI;
	get_rotation_matrix(axis, t, orient);
}

void matrix_matrix_multiplication(matrix m, matrix n, matrix res) {
	res[0][0] = m[0][0]*n[0][0] + m[0][1]*n[1][0] + m[0][2]*n[2][0];
	res[0][1] = m[0][0]*n[0][1] + m[0][1]*n[1][1] + m[0][2]*n[2][1];
	res[0][2] = m[0][0]*n[0][2] + m[0][1]*n[1][2] + m[0][2]*n[2][2];

	res[1][0] = m[1][0]*n[0][0] + m[1][1]*n[1][0] + m[1][2]*n[2][0];
	res[1][1] = m[1][0]*n[0][1] + m[1][1]*n[1][1] + m[1][2]*n[2][1];
	res[1][2] = m[1][0]*n[0][2] + m[1][1]*n[1][2] + m[1][2]*n[2][2];

	res[2][0] = m[2][0]*n[0][0] + m[2][1]*n[1][0] + m[2][2]*n[2][0];
	res[2][1] = m[2][0]*n[0][1] + m[2][1]*n[1][1] + m[2][2]*n[2][1];
	res[2][2] = m[2][0]*n[0][2] + m[2][1]*n[1][2] + m[2][2]*n[2][2];
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

void utils_rotate_matrix(matrix orient_old, matrix orient_new, vector axis, double t) {
	matrix rotation_matrix;
	get_rotation_matrix(axis, t, rotation_matrix);

	matrix_matrix_multiplication(orient_old, rotation_matrix, orient_new);
}

void get_rotated_vector(vector v, vector axis, double t, vector res) {
	matrix rotation_matrix;
	get_rotation_matrix(axis, t, rotation_matrix);
	MATRIX_VECTOR_MULTIPLICATION(rotation_matrix, v, res);
}

void rotate_vector(vector v, vector axis, double t) {
	vector tmp;
	get_rotated_vector(v, axis, t, tmp);
	v[0] = tmp[0];
	v[1] = tmp[1];
	v[2] = tmp[2];
}

void set_orientation_around_vector(vector v, matrix orient, double t) {
	vector y_tmp, w_tmp, z;

	y_tmp[0] = -v[0];
	y_tmp[1] = -v[1];
	y_tmp[2] = -v[2];

	random_vector_on_sphere(z);
	random_vector_on_sphere(w_tmp);

	// orthonormalize
	gram_schmidt(y_tmp, z, w_tmp);

	matrix rotation_matrix;
	get_rotation_matrix(z, t, rotation_matrix);

	vector y, w;
	MATRIX_VECTOR_MULTIPLICATION(rotation_matrix, y_tmp, y);
	MATRIX_VECTOR_MULTIPLICATION(rotation_matrix, w_tmp, w);

	matrix cambiamento_base;

	cambiamento_base[0][0] = z[0];
	cambiamento_base[1][0] = z[1];
	cambiamento_base[2][0] = z[2];

	cambiamento_base[0][1] = w[0];
	cambiamento_base[1][1] = w[1];
	cambiamento_base[2][1] = w[2];

	cambiamento_base[0][2] = y[0];
	cambiamento_base[1][2] = y[1];
	cambiamento_base[2][2] = y[2];

	// rotations have det(R) == 1
	if(determinant(cambiamento_base) < 0) {
		cambiamento_base[0][0] = w[0];
		cambiamento_base[1][0] = w[1];
		cambiamento_base[2][0] = w[2];

		cambiamento_base[0][1] = z[0];
		cambiamento_base[1][1] = z[1];
		cambiamento_base[2][1] = z[2];
	}

	int i, j;
	matrix orient_old;
	for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) orient_old[i][j] = orient[i][j];
	matrix_matrix_multiplication(cambiamento_base, orient_old, orient);
}

void utils_reset_acceptance_counters(System *syst) {
	int i;
	for(i = 0; i < N_MOVES; i++) {
		syst->tries[i] = 0;
		syst->accepted[i] = 0;
	}
}


int getLine(char *line,FILE *pfile)
{
	if (fgets(line,MAX_LINE_LENGTH,pfile) == NULL)
		return 0;
	else
		return strlen(line);
}
