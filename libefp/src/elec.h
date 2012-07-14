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

#ifndef LIBEFP_ELEC_H
#define LIBEFP_ELEC_H

static inline int
quad_idx(int a, int b)
{
	/* order in which GAMESS stores quadrupoles */
	enum { xx = 0, yy, zz, xy, xz, yz };

	static const int idx[] = {
		xx, xy, xz, xy, yy, yz, xz, yz, zz
	};

	return idx[a * 3 + b];
}

static inline int
oct_idx(int a, int b, int c)
{
	/* order in which GAMESS stores octupoles */
	enum { xxx = 0, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz };

	static const int idx[] = {
		xxx, xxy, xxz, xxy, xyy, xyz, xxz, xyz, xzz,
		xxy, xyy, xyz, xyy, yyy, yyz, xyz, yyz, yzz,
		xxz, xyz, xzz, xyz, yyz, yzz, xzz, yzz, zzz
	};

	return idx[a * 9 + b * 3 + c];
}

static inline double
quadrupole_sum(const double *quad, const vec_t *dr)
{
	/* order in which GAMESS stores quadrupoles */
	enum { xx = 0, yy, zz, xy, xz, yz };

	double sum = 0.0;

	sum += quad[xx] * dr->x * dr->x;
	sum += quad[yy] * dr->y * dr->y;
	sum += quad[zz] * dr->z * dr->z;
	sum += quad[xy] * dr->x * dr->y * 2.0;
	sum += quad[xz] * dr->x * dr->z * 2.0;
	sum += quad[yz] * dr->y * dr->z * 2.0;

	return sum;
}

static inline double
octupole_sum(const double *oct, const vec_t *dr)
{
	/* order in which GAMESS stores octupoles */
	enum { xxx = 0, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz };

	double sum = 0.0;

	sum += oct[xxx] * dr->x * dr->x * dr->x;
	sum += oct[yyy] * dr->y * dr->y * dr->y;
	sum += oct[zzz] * dr->z * dr->z * dr->z;
	sum += oct[xxy] * dr->x * dr->x * dr->y * 3.0;
	sum += oct[xxz] * dr->x * dr->x * dr->z * 3.0;
	sum += oct[xyy] * dr->x * dr->y * dr->y * 3.0;
	sum += oct[yyz] * dr->y * dr->y * dr->z * 3.0;
	sum += oct[xzz] * dr->x * dr->z * dr->z * 3.0;
	sum += oct[yzz] * dr->y * dr->z * dr->z * 3.0;
	sum += oct[xyz] * dr->x * dr->y * dr->z * 6.0;

	return sum;
}

void efp_charge_charge_grad(double q1, double q2, const vec_t *dr,
			    vec_t *force, vec_t *add1, vec_t *add2);

void efp_charge_dipole_grad(double q1, const vec_t *d2, const vec_t *dr,
			    vec_t *force, vec_t *add1, vec_t *add2);

void efp_charge_quadrupole_grad(double q1, const double *quad2, const vec_t *dr,
				vec_t *force, vec_t *add1, vec_t *add2);

void efp_charge_octupole_grad(double q1, const double *oct2, const vec_t *dr,
			      vec_t *force, vec_t *add1, vec_t *add2);

void efp_dipole_dipole_grad(const vec_t *d1, const vec_t *d2, const vec_t *dr,
			    vec_t *force, vec_t *add1, vec_t *add2);

void efp_dipole_quadrupole_grad(const vec_t *d1, const double *quad2,
				const vec_t *dr, vec_t *force, vec_t *add1,
				vec_t *add2);

void efp_quadrupole_quadrupole_grad(const double *quad1, const double *quad2,
				    const vec_t *dr, vec_t *force, vec_t *add1,
				    vec_t *add2);

#endif /* LIBEFP_ELEC_H */
