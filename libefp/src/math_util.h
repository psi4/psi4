/*-
 * Copyright (c) 2012 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef LIBEFP_MATH_UTIL_H
#define LIBEFP_MATH_UTIL_H

#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define VEC(x) ((vec_t *)(&(x)))

typedef struct {
	double x, y, z;
} vec_t;

typedef struct {
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
} mat_t;

static inline double
vec_el(const vec_t *vec, int idx)
{
	return ((double *)vec)[idx];
}

static inline void
vec_zero(vec_t *vec)
{
	vec->x = 0.0, vec->y = 0.0, vec->z = 0.0;
}

static inline void
vec_negate(vec_t *vec)
{
	vec->x = -vec->x, vec->y = -vec->y, vec->z = -vec->z;
}

static inline double
vec_dot(const vec_t *a, const vec_t *b)
{
	return a->x * b->x + a->y * b->y + a->z * b->z;
}

static inline vec_t
vec_cross(const vec_t *a, const vec_t *b)
{
	vec_t c = {
		a->y * b->z - a->z * b->y,
		a->z * b->x - a->x * b->z,
		a->x * b->y - a->y * b->x
	};
	return c;
}

static inline vec_t
vec_add(const vec_t *a, const vec_t *b)
{
	vec_t c = { a->x + b->x, a->y + b->y, a->z + b->z };
	return c;
}

static inline vec_t
vec_sub(const vec_t *a, const vec_t *b)
{
	vec_t c = { a->x - b->x, a->y - b->y, a->z - b->z };
	return c;
}

static inline double
vec_len_2(const vec_t *a)
{
	return vec_dot(a, a);
}

static inline double
vec_len(const vec_t *a)
{
	return sqrt(vec_len_2(a));
}

static inline double
vec_dist_2(const vec_t *a, const vec_t *b)
{
	vec_t dr = vec_sub(a, b);
	return vec_len_2(&dr);
}

static inline double
vec_dist(const vec_t *a, const vec_t *b)
{
	return sqrt(vec_dist_2(a, b));
}

static inline int
eq(double a, double b)
{
	static const double eps = 1.0e-8;
	return fabs(a - b) < eps;
}

static inline void
mat_zero(mat_t *mat)
{
	mat->xx = 0.0, mat->xy = 0.0, mat->xz = 0.0;
	mat->yx = 0.0, mat->yy = 0.0, mat->yz = 0.0;
	mat->zx = 0.0, mat->zy = 0.0, mat->zz = 0.0;
}

static inline void
mat_vec(const mat_t *mat, const vec_t *vec, vec_t *out)
{
	out->x = mat->xx * vec->x + mat->xy * vec->y + mat->xz * vec->z;
	out->y = mat->yx * vec->x + mat->yy * vec->y + mat->yz * vec->z;
	out->z = mat->zx * vec->x + mat->zy * vec->y + mat->zz * vec->z;
}

static inline void
mat_trans_vec(const mat_t *mat, const vec_t *vec, vec_t *out)
{
	out->x = mat->xx * vec->x + mat->yx * vec->y + mat->zx * vec->z;
	out->y = mat->xy * vec->x + mat->yy * vec->y + mat->zy * vec->z;
	out->z = mat->xz * vec->x + mat->yz * vec->y + mat->zz * vec->z;
}

static inline void
move_pt(const vec_t *com, const mat_t *rotmat,
	const vec_t *com_init, const vec_t *pos_init,
	vec_t *out)
{
	vec_t pos = vec_sub(pos_init, com_init);
	mat_vec(rotmat, &pos, out);
	out->x += com->x, out->y += com->y, out->z += com->z;
}

static inline void
rotate_t2(const mat_t *rotmat, const double *in, double *out)
{
	const double *rm = (const double *)rotmat;
	memset(out, 0, 9 * sizeof(double));

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
			out[a2 * 3 + b2] += in[a1 * 3 + b1] *
					rm[a2 * 3 + a1] *
					rm[b2 * 3 + b1];
}

static inline void
rotate_t3(const mat_t *rotmat, const double *in, double *out)
{
	const double *rm = (const double *)rotmat;
	memset(out, 0, 27 * sizeof(double));

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
	for (int c1 = 0; c1 < 3; c1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
		for (int c2 = 0; c2 < 3; c2++)
			out[a2 * 9 + b2 * 3 + c2] += in[a1 * 9 + b1 * 3 + c1] *
					rm[a2 * 3 + a1] *
					rm[b2 * 3 + b1] *
					rm[c2 * 3 + c1];
}

#endif /* LIBEFP_MATH_UTIL_H */
