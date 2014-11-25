/*-
 * Copyright (c) 2012-2014 Ilya Kaliman
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

#include "math_util.h"

static inline void
add_3(vec_t *a, const vec_t *aa,
      vec_t *b, const vec_t *bb,
      vec_t *c, const vec_t *cc)
{
	a->x += aa->x;
	a->y += aa->y;
	a->z += aa->z;

	b->x += bb->x;
	b->y += bb->y;
	b->z += bb->z;

	c->x += cc->x;
	c->y += cc->y;
	c->z += cc->z;
}

static inline size_t
quad_idx(size_t a, size_t b)
{
	/* order in which quadrupoles are stored */
	enum { xx = 0, yy, zz, xy, xz, yz };

	static const size_t idx[] = {
		xx, xy, xz, xy, yy, yz, xz, yz, zz
	};

	return idx[a * 3 + b];
}

static inline size_t
oct_idx(size_t a, size_t b, size_t c)
{
	/* order in which octupoles are stored */
	enum { xxx = 0, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz };

	static const size_t idx[] = {
		xxx, xxy, xxz, xxy, xyy, xyz, xxz, xyz, xzz,
		xxy, xyy, xyz, xyy, yyy, yyz, xyz, yyz, yzz,
		xxz, xyz, xzz, xyz, yyz, yzz, xzz, yzz, zzz
	};

	return idx[a * 9 + b * 3 + c];
}

static inline double
quadrupole_sum(const double *quad, const vec_t *dr)
{
	/* order in which quadrupoles are stored */
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

double efp_charge_charge_energy(double, double, const vec_t *);
double efp_charge_dipole_energy(double, const vec_t *, const vec_t *);
double efp_charge_quadrupole_energy(double, const double *, const vec_t *);
double efp_charge_octupole_energy(double, const double *, const vec_t *);
double efp_dipole_dipole_energy(const vec_t *, const vec_t *, const vec_t *);
double efp_dipole_quadrupole_energy(const vec_t *, const double *, const vec_t *);
double efp_quadrupole_quadrupole_energy(const double *, const double *, const vec_t *);

void efp_charge_charge_grad(double, double, const vec_t *, vec_t *, vec_t *, vec_t *);
void efp_charge_dipole_grad(double, const vec_t *, const vec_t *, vec_t *, vec_t *, vec_t *);
void efp_charge_quadrupole_grad(double, const double *, const vec_t *, vec_t *, vec_t *, vec_t *);
void efp_charge_octupole_grad(double, const double *, const vec_t *, vec_t *, vec_t *, vec_t *);
void efp_dipole_dipole_grad(const vec_t *, const vec_t *, const vec_t *, vec_t *, vec_t *, vec_t *);
void efp_dipole_quadrupole_grad(const vec_t *, const double *, const vec_t *, vec_t *, vec_t *, vec_t *);
void efp_quadrupole_quadrupole_grad(const double *, const double *, const vec_t *, vec_t *, vec_t *, vec_t *);

#endif /* LIBEFP_ELEC_H */
