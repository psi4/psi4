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

#include "math_util.h"
#include "elec.h"

static double
octupole_sum_xyz(const double *oct, const vec_t *dr, int axis)
{
	const double *pdr = (const double *)dr;
	double sum = 0.0;

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				double o = oct[oct_idx(a, b, c)];
				if (a == axis) sum += o * pdr[b] * pdr[c];
				if (b == axis) sum += o * pdr[a] * pdr[c];
				if (c == axis) sum += o * pdr[a] * pdr[b];
			}

	return sum;
}

void
efp_charge_charge_grad(double q1, double q2, const vec_t *dr,
		       vec_t *force, vec_t *add1, vec_t *add2)
{
	double r = vec_len(dr);
	double r3 = r * r * r;

	double g = q1 * q2 / r3;

	force->x = g * dr->x;
	force->y = g * dr->y;
	force->z = g * dr->z;

	add1->x = 0.0;
	add1->y = 0.0;
	add1->z = 0.0;

	add2->x = 0.0;
	add2->y = 0.0;
	add2->z = 0.0;
}

void
efp_charge_dipole_grad(double q1, const vec_t *d2, const vec_t *dr,
		       vec_t *force, vec_t *add1, vec_t *add2)
{
	double r = vec_len(dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;

	double t1 = 3.0 * q1 / r5 * vec_dot(d2, dr);
	double t2 = q1 / r3;

	force->x = t2 * d2->x - t1 * dr->x;
	force->y = t2 * d2->y - t1 * dr->y;
	force->z = t2 * d2->z - t1 * dr->z;

	add1->x = 0.0;
	add1->y = 0.0;
	add1->z = 0.0;

	add2->x = t2 * (d2->y * dr->z - d2->z * dr->y);
	add2->y = t2 * (d2->z * dr->x - d2->x * dr->z);
	add2->z = t2 * (d2->x * dr->y - d2->y * dr->x);
}

void
efp_charge_quadrupole_grad(double q1, const double *quad2, const vec_t *dr,
			   vec_t *force, vec_t *add1, vec_t *add2)
{
	double r = vec_len(dr);
	double r2 = r * r;
	double r5 = r2 * r2 * r;
	double r7 = r5 * r2;

	double t1x = q1 / r5 * -2.0 * (dr->x * quad2[quad_idx(0, 0)] +
				       dr->y * quad2[quad_idx(0, 1)] +
				       dr->z * quad2[quad_idx(0, 2)]);
	double t1y = q1 / r5 * -2.0 * (dr->x * quad2[quad_idx(1, 0)] +
				       dr->y * quad2[quad_idx(1, 1)] +
				       dr->z * quad2[quad_idx(1, 2)]);
	double t1z = q1 / r5 * -2.0 * (dr->x * quad2[quad_idx(2, 0)] +
				       dr->y * quad2[quad_idx(2, 1)] +
				       dr->z * quad2[quad_idx(2, 2)]);

	double g = 5.0 * q1 / r7 * quadrupole_sum(quad2, dr);

	force->x = g * dr->x + t1x;
	force->y = g * dr->y + t1y;
	force->z = g * dr->z + t1z;

	add1->x = 0.0;
	add1->y = 0.0;
	add1->z = 0.0;

	add2->x = t1z * dr->y - t1y * dr->z;
	add2->y = t1x * dr->z - t1z * dr->x;
	add2->z = t1y * dr->x - t1x * dr->y;
}

void
efp_charge_octupole_grad(double q1, const double *oct2, const vec_t *dr,
			 vec_t *force, vec_t *add1, vec_t *add2)
{
	double r = vec_len(dr);
	double r3 = r * r * r;
	double r7 = r3 * r3 * r;
	double r9 = r3 * r3 * r3;

	double t1x = q1 / r7 * octupole_sum_xyz(oct2, dr, 0);
	double t1y = q1 / r7 * octupole_sum_xyz(oct2, dr, 1);
	double t1z = q1 / r7 * octupole_sum_xyz(oct2, dr, 2);

	double g = 7.0 * q1 / r9 * octupole_sum(oct2, dr);

	force->x = -g * dr->x + t1x;
	force->y = -g * dr->y + t1y;
	force->z = -g * dr->z + t1z;

	add1->x = 0.0;
	add1->y = 0.0;
	add1->z = 0.0;

	add2->x = t1y * dr->z - t1z * dr->y;
	add2->y = t1z * dr->x - t1x * dr->z;
	add2->z = t1x * dr->y - t1y * dr->x;
}

void
efp_dipole_dipole_grad(const vec_t *d1, const vec_t *d2, const vec_t *dr,
		       vec_t *force, vec_t *add1, vec_t *add2)
{
	double r = vec_len(dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;
	double r7 = r5 * r * r;

	double d1dr = vec_dot(d1, dr);
	double d2dr = vec_dot(d2, dr);

	double t1 = 3.0 / r5;
	double t2 = t1 * vec_dot(d1, d2) - 15.0 / r7 * d1dr * d2dr;

	force->x = t2 * dr->x + t1 * (d2dr * d1->x + d1dr * d2->x);
	force->y = t2 * dr->y + t1 * (d2dr * d1->y + d1dr * d2->y);
	force->z = t2 * dr->z + t1 * (d2dr * d1->z + d1dr * d2->z);

	add1->x = d1->y * (d2->z / r3 - t1 * dr->z * d2dr) -
			d1->z * (d2->y / r3 - t1 * dr->y * d2dr);
	add1->y = d1->z * (d2->x / r3 - t1 * dr->x * d2dr) -
			d1->x * (d2->z / r3 - t1 * dr->z * d2dr);
	add1->z = d1->x * (d2->y / r3 - t1 * dr->y * d2dr) -
			d1->y * (d2->x / r3 - t1 * dr->x * d2dr);

	add2->x = d2->y * (d1->z / r3 - t1 * dr->z * d1dr) -
			d2->z * (d1->y / r3 - t1 * dr->y * d1dr);
	add2->y = d2->z * (d1->x / r3 - t1 * dr->x * d1dr) -
			d2->x * (d1->z / r3 - t1 * dr->z * d1dr);
	add2->z = d2->x * (d1->y / r3 - t1 * dr->y * d1dr) -
			d2->y * (d1->x / r3 - t1 * dr->x * d1dr);
}

void
efp_dipole_quadrupole_grad(const vec_t *d1, const double *quad2,
			   const vec_t *dr, vec_t *force, vec_t *add1,
			   vec_t *add2)
{
	double r = vec_len(dr);
	double r2 = r * r;
	double r3 = r2 * r;
	double r5 = r3 * r2;
	double r7 = r5 * r2;
	double r9 = r7 * r2;

	double q2sx = 0.0;
	double q2sy = 0.0;
	double q2sz = 0.0;

	for (int a = 0; a < 3; a++) {
		q2sx += quad2[quad_idx(0, a)] * vec_el(dr, a);
		q2sy += quad2[quad_idx(1, a)] * vec_el(dr, a);
		q2sz += quad2[quad_idx(2, a)] * vec_el(dr, a);
	}

	double d1dr = vec_dot(d1, dr);
	double q2s = quadrupole_sum(quad2, dr);

	double t1 = d1->x * q2sx + d1->y * q2sy + d1->z * q2sz;
	double t2 = -10.0 / r7 * t1 + 35.0 / r9 * q2s * d1dr;

	double d1q2x = d1->x * quad2[quad_idx(0, 0)] +
		       d1->y * quad2[quad_idx(0, 1)] +
		       d1->z * quad2[quad_idx(0, 2)];
	double d1q2y = d1->x * quad2[quad_idx(1, 0)] +
		       d1->y * quad2[quad_idx(1, 1)] +
		       d1->z * quad2[quad_idx(1, 2)];
	double d1q2z = d1->x * quad2[quad_idx(2, 0)] +
		       d1->y * quad2[quad_idx(2, 1)] +
		       d1->z * quad2[quad_idx(2, 2)];

	double q2xdr = dr->x * quad2[quad_idx(0, 0)] +
		       dr->y * quad2[quad_idx(0, 1)] +
		       dr->z * quad2[quad_idx(0, 2)];
	double q2ydr = dr->x * quad2[quad_idx(1, 0)] +
		       dr->y * quad2[quad_idx(1, 1)] +
		       dr->z * quad2[quad_idx(1, 2)];
	double q2zdr = dr->x * quad2[quad_idx(2, 0)] +
		       dr->y * quad2[quad_idx(2, 1)] +
		       dr->z * quad2[quad_idx(2, 2)];

	force->x = t2 * dr->x + 2.0 / r5 * d1q2x -
			5.0 / r7 * (q2s * d1->x + 2.0 * q2xdr * d1dr);
	force->y = t2 * dr->y + 2.0 / r5 * d1q2y -
			5.0 / r7 * (q2s * d1->y + 2.0 * q2ydr * d1dr);
	force->z = t2 * dr->z + 2.0 / r5 * d1q2z -
			5.0 / r7 * (q2s * d1->z + 2.0 * q2zdr * d1dr);

	add1->x = 2.0 / r5 * (d1->z * q2ydr - d1->y * q2zdr) +
			5.0 / r7 * q2s * (dr->z * d1->y - dr->y * d1->z);
	add1->y = 2.0 / r5 * (d1->x * q2zdr - d1->z * q2xdr) +
			5.0 / r7 * q2s * (dr->x * d1->z - dr->z * d1->x);
	add1->z = 2.0 / r5 * (d1->y * q2xdr - d1->x * q2ydr) +
			5.0 / r7 * q2s * (dr->y * d1->x - dr->x * d1->y);

	add2->x = -10.0 / r7 * d1dr * (q2ydr * dr->z - q2zdr * dr->y) -
			2.0 / r5 * ((q2zdr * d1->y + dr->y * d1q2z) -
			      (q2ydr * d1->z + dr->z * d1q2y));
	add2->y = -10.0 / r7 * d1dr * (q2zdr * dr->x - q2xdr * dr->z) -
			2.0 / r5 * ((q2xdr * d1->z + dr->z * d1q2x) -
			      (q2zdr * d1->x + dr->x * d1q2z));
	add2->z = -10.0 / r7 * d1dr * (q2xdr * dr->y - q2ydr * dr->x) -
			2.0 / r5 * ((q2ydr * d1->x + dr->x * d1q2y) -
			      (q2xdr * d1->y + dr->y * d1q2x));
}

void
efp_quadrupole_quadrupole_grad(const double *quad1, const double *quad2,
			       const vec_t *dr, vec_t *force, vec_t *add1,
			       vec_t *add2)
{
	double r = vec_len(dr);
	double r2 = r * r;
	double r5 = r2 * r2 * r;
	double r7 = r5 * r2;
	double r9 = r7 * r2;
	double r11 = r9 * r2;

	double q1ss = quadrupole_sum(quad1, dr);
	double q2ss = quadrupole_sum(quad2, dr);

	double q1s[3] = { 0.0, 0.0, 0.0 };
	double q2s[3] = { 0.0, 0.0, 0.0 };

	double q1sq2s = 0.0;

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			q1s[a] += quad1[quad_idx(a, b)] * vec_el(dr, b);
			q2s[a] += quad2[quad_idx(a, b)] * vec_el(dr, b);
		}
		q1sq2s += q1s[a] * q2s[a];
	}

	double q1q2 = 0.0;

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			q1q2 += quad1[quad_idx(a, b)] * quad2[quad_idx(a, b)];

	double g = 30.0 / r7 * q1q2 - 420.0 / r9 * q1sq2s +
			945.0 / r11 * q1ss * q2ss;

	double t1x = 0.0, t1y = 0.0, t1z = 0.0;

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++) {
			int ab = quad_idx(a, b);
			double dra = vec_el(dr, a);
			t1x += (quad1[quad_idx(0, b)] * quad2[ab] +
				quad1[ab] * quad2[quad_idx(0, b)]) * dra;
			t1y += (quad1[quad_idx(1, b)] * quad2[ab] +
				quad1[ab] * quad2[quad_idx(1, b)]) * dra;
			t1z += (quad1[quad_idx(2, b)] * quad2[ab] +
				quad1[ab] * quad2[quad_idx(2, b)]) * dra;
		}

	force->x = (g * dr->x + 60.0 / r7 * t1x -
			210.0 / r9 * (q1s[0] * q2ss + q2s[0] * q1ss)) / 9.0;
	force->y = (g * dr->y + 60.0 / r7 * t1y -
			210.0 / r9 * (q1s[1] * q2ss + q2s[1] * q1ss)) / 9.0;
	force->z = (g * dr->z + 60.0 / r7 * t1z -
			210.0 / r9 * (q1s[2] * q2ss + q2s[2] * q1ss)) / 9.0;

	double q1q2tt[3][3];
	memset(q1q2tt, 0, 9 * sizeof(double));

	for (int a = 0; a < 3; a++)
	for (int b = 0; b < 3; b++)
	for (int c = 0; c < 3; c++) {
		double dra = vec_el(dr, a);
		double drc = vec_el(dr, c);
		q1q2tt[b][c] += quad1[quad_idx(a, b)] *
				(-10.0 / r7 * (drc * q2s[a] + dra * q2s[c]) +
				  35.0 / r9 * dra * drc * q2ss +
				   2.0 / r5 * quad2[quad_idx(a, c)]);
	}

	add1->x = 2.0 / 3.0 * (q1q2tt[1][2] - q1q2tt[2][1]);
	add1->y = 2.0 / 3.0 * (q1q2tt[2][0] - q1q2tt[0][2]);
	add1->z = 2.0 / 3.0 * (q1q2tt[0][1] - q1q2tt[1][0]);

	double q2q1tt[3][3];
	memset(q2q1tt, 0, 9 * sizeof(double));

	for (int a = 0; a < 3; a++)
	for (int b = 0; b < 3; b++)
	for (int c = 0; c < 3; c++) {
		double dra = vec_el(dr, a);
		double drc = vec_el(dr, c);
		q2q1tt[b][c] += quad2[quad_idx(a, b)] *
				(-10.0 / r7 * (drc * q1s[a] + dra * q1s[c]) +
				  35.0 / r9 * dra * drc * q1ss +
				   2.0 / r5 * quad1[quad_idx(a, c)]);
	}

	add2->x = 2.0 / 3.0 * (q2q1tt[1][2] - q2q1tt[2][1]);
	add2->y = 2.0 / 3.0 * (q2q1tt[2][0] - q2q1tt[0][2]);
	add2->z = 2.0 / 3.0 * (q2q1tt[0][1] - q2q1tt[1][0]);
}
