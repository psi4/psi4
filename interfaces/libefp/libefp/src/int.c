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

#include <assert.h>
#include <string.h>

#include "int.h"
#include "int_shift.h"

/* Overlap and kinetic energy integral computation routines. */

/* Tolerance for integrals (20 * ln10) */
static const double int_tol = 46.051701859881;

/* Normalization constants */
static const double int_norm[] = {
	1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 1.0000000000000000,
	1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 1.7320508075688801,
	1.7320508075688801, 1.7320508075688801, 1.0000000000000000, 1.0000000000000000,
	1.0000000000000000, 2.2360679774997898, 2.2360679774997898, 2.2360679774997898,
	2.2360679774997898, 2.2360679774997898, 2.2360679774997898, 3.8729833462074232
};

static void set_coef(double *con, char type, const double *coef)
{
	switch (type) {
		case 'S':
			con[0] = *coef;
			return;
		case 'L':
			con[0] = *coef++;
			/* fall through */
		case 'P':
			for (size_t i = 1; i < 4; i++)
				con[i] = *coef;
			return;
		case 'D':
			for (size_t i = 4; i < 10; i++)
				con[i] = *coef;
			return;
		case 'F':
			for (size_t i = 10; i < 20; i++)
				con[i] = *coef;
			return;
	}
}

static void make_int(size_t ni, size_t nj, double tt, const vec_t *p,
		     const vec_t *p_i, const vec_t *p_j, vec_t *out)
{
	static const size_t imin[] = { 0, 1, 3,  6, 10, 15, 21, 28, 36, 45 };
	static const size_t imax[] = { 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 };

	static const double h[] = {
		 0.0000000000000000, -0.7071067811865475,  0.7071067811865475,
		-1.2247448713915889,  0.0000000000000000,  1.2247448713915889,
		-1.6506801238857844, -0.5246476232752903,  0.5246476232752903,
		 1.6506801238857844, -2.0201828704560856, -0.9585724646138185,
		 0.0000000000000000,  0.9585724646138185,  2.0201828704560856,
		-2.3506049736744923, -1.3358490740136970, -0.4360774119276165,
		 0.4360774119276165,  1.3358490740136970,  2.3506049736744923,
		-2.6519613568352334, -1.6735516287674714, -0.8162878828589647,
		 0.0000000000000000,  0.8162878828589647,  1.6735516287674714,
		 2.6519613568352334
	};

	static const double w[] = {
		 1.7724538509055161,  0.8862269254527580,  0.8862269254527580,
		 0.2954089751509193,  1.1816359006036774,  0.2954089751509193,
		 0.0813128354472451,  0.8049140900055128,  0.8049140900055128,
		 0.0813128354472451,  0.0199532420590459,  0.3936193231522411,
		 0.9453087204829419,  0.3936193231522411,  0.0199532420590459,
		 0.0045300099055088,  0.1570673203228566,  0.7246295952243925,
		 0.7246295952243925,  0.1570673203228566,  0.0045300099055088,
		 0.0009717812450995,  0.0545155828191270,  0.4256072526101277,
		 0.8102646175568073,  0.4256072526101277,  0.0545155828191270,
		 0.0009717812450995
	};

	size_t npts = (ni + nj) / 2;

	double xint = 0.0, yint = 0.0, zint = 0.0;

	for (size_t i = imin[npts]; i < imax[npts]; i++) {
		double px = w[i];
		double py = w[i];
		double pz = w[i];

		double tmp = h[i] * tt;

		if(ni > 0) {
			double ax = tmp + p->x - p_i->x;
			double ay = tmp + p->y - p_i->y;
			double az = tmp + p->z - p_i->z;

			/* fancy loop unrolling */
			switch (ni) {
				case 4:
					px *= ax;
					py *= ay;
					pz *= az;
				case 3:
					px *= ax;
					py *= ay;
					pz *= az;
				case 2:
					px *= ax;
					py *= ay;
					pz *= az;
				case 1:
					px *= ax;
					py *= ay;
					pz *= az;
					break;
				default:
					assert(0);
			}
		}

		if(nj > 0) {
			double bx = tmp + p->x - p_j->x;
			double by = tmp + p->y - p_j->y;
			double bz = tmp + p->z - p_j->z;

			/* fancy loop unrolling */
			switch (nj) {
				case 5:
					px *= bx;
					py *= by;
					pz *= bz;
				case 4:
					px *= bx;
					py *= by;
					pz *= bz;
				case 3:
					px *= bx;
					py *= by;
					pz *= bz;
				case 2:
					px *= bx;
					py *= by;
					pz *= bz;
				case 1:
					px *= bx;
					py *= by;
					pz *= bz;
					break;
				default:
					assert(0);
			}
		}

		xint += px;
		yint += py;
		zint += pz;
	}

	out->x = xint;
	out->y = yint;
	out->z = zint;
}

static size_t get_shell_idx(char type)
{
	switch (type) {
		case 'S': return 0;
		case 'L': return 1;
		case 'P': return 2;
		case 'D': return 3;
		case 'F': return 4;
	}
	assert(0);
}

static size_t get_shell_start(size_t shell_idx)
{
	static const size_t start[] = {
		0, 0, 1, 4, 10
	};

	return start[shell_idx];
}

static size_t get_shell_end(size_t shell_idx)
{
	static const size_t end[] = {
		1, 4, 4, 10, 20
	};

	return end[shell_idx];
}

static size_t get_shell_sl(size_t shell_idx)
{
	static const size_t sl[] = {
		1, 2, 2, 3, 4
	};

	return sl[shell_idx];
}

static void init_ft(size_t count_i, char type_j, double *ft)
{
	switch (type_j) {
		case 'S':
			for (size_t i = 0; i < count_i; i++) {
				*ft++ = 3.0;
			}
			return;
		case 'L':
			for (size_t i = 0; i < count_i; i++) {
				*ft++ = 3.0;
				*ft++ = 5.0;
				*ft++ = 5.0;
				*ft++ = 5.0;
			}
			return;
		case 'P':
			for (size_t i = 0; i < count_i; i++) {
				*ft++ = 5.0;
				*ft++ = 5.0;
				*ft++ = 5.0;
			}
			return;
		case 'D':
			for (size_t i = 0; i < count_i; i++) {
				*ft++ = 7.0;
				*ft++ = 7.0;
				*ft++ = 7.0;
				*ft++ = 7.0;
				*ft++ = 7.0;
				*ft++ = 7.0;
			}
			return;
		case 'F':
			for (size_t i = 0; i < count_i; i++) {
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
				*ft++ = 9.0;
			}
			return;
	}
	assert(0);
}

void efp_st_int(size_t n_atoms_i, const struct xr_atom *atoms_i,
		size_t n_atoms_j, const struct xr_atom *atoms_j,
		size_t stride, double *s, double *t)
{
	double xin[90];
	double yin[90];
	double zin[90];

	double ft[100], dij[100];
	double sblk[100], tblk[100];

	for (size_t iii = 0, loc_i = 0; iii < n_atoms_i; iii++) {
		const struct xr_atom *at_i = atoms_i + iii;

	/* shell i */
	for (size_t ii = 0; ii < at_i->n_shells; ii++) {
		const struct shell *sh_i = at_i->shells + ii;

		size_t type_i = get_shell_idx(sh_i->type);
		size_t start_i = get_shell_start(type_i);
		size_t end_i = get_shell_end(type_i);
		size_t sl_i = get_shell_sl(type_i);
		size_t count_i = end_i - start_i;

		for (size_t jjj = 0, loc_j = 0; jjj < n_atoms_j; jjj++) {
			const struct xr_atom *at_j = atoms_j + jjj;

		/* shell j */
		for (size_t jj = 0; jj < at_j->n_shells; jj++) {
			const struct shell *sh_j = at_j->shells + jj;

			size_t type_j = get_shell_idx(sh_j->type);
			size_t start_j = get_shell_start(type_j);
			size_t end_j = get_shell_end(type_j);
			size_t sl_j = get_shell_sl(type_j);
			size_t count_j = end_j - start_j;

			size_t count = count_i * count_j;

			memset(sblk, 0, count * sizeof(double));
			memset(tblk, 0, count * sizeof(double));

			init_ft(count_i, sh_j->type, ft);

			double rr = vec_dist_2(CVEC(at_i->x), CVEC(at_j->x));

			const size_t *shift_x = shift_table_x[type_i * 5 + type_j];
			const size_t *shift_y = shift_table_y[type_i * 5 + type_j];
			const size_t *shift_z = shift_table_z[type_i * 5 + type_j];

			const double *coef_i = sh_i->coef;

			/* primitive i */
			for (size_t ig = 0; ig < sh_i->n_funcs; ig++) {
				double ai = *coef_i++;

				double con_i[20];
				set_coef(con_i, sh_i->type, coef_i);

				coef_i++;
				if (sh_i->type == 'L')
					coef_i++;

				const double *coef_j = sh_j->coef;

				/* primitive j */
				for (size_t jg = 0; jg < sh_j->n_funcs; jg++) {
					double aj = *coef_j++;

					double aa = 1.0 / (ai + aj);
					double tmp = aj * ai * rr * aa;

					if (tmp > int_tol) {
						coef_j++;
						if (sh_j->type == 'L')
							coef_j++;

						continue;
					}

					double con_j[20];
					set_coef(con_j, sh_j->type, coef_j);

					coef_j++;
					if (sh_j->type == 'L')
						coef_j++;

					vec_t a = {
						(ai * at_i->x + aj * at_j->x) * aa,
						(ai * at_i->y + aj * at_j->y) * aa,
						(ai * at_i->z + aj * at_j->z) * aa
					};

					double fac = exp(-tmp);

					for (size_t i = start_i, idx = 0; i < end_i; i++)
						for (size_t j = start_j; j < end_j; j++, idx++)
							dij[idx] = fac * con_i[i] * int_norm[i] * con_j[j] * int_norm[j];

					double taa = sqrt(aa);
					double t1 = -2.0 * aj * aj * taa;
					double t2 = -0.5 * taa;

					for (size_t i = 0, idx = 0; i < sl_i; i++, idx += 5) {
						for (size_t j = 0; j < sl_j; j++) {
							vec_t iout;

							make_int(i, j, taa, &a, CVEC(at_i->x), CVEC(at_j->x), &iout);
							xin[idx + j] = iout.x * taa;
							yin[idx + j] = iout.y * taa;
							zin[idx + j] = iout.z * taa;

							make_int(i, j + 2, taa, &a, CVEC(at_i->x), CVEC(at_j->x), &iout);
							xin[idx + j + 30] = iout.x * t1;
							yin[idx + j + 30] = iout.y * t1;
							zin[idx + j + 30] = iout.z * t1;

							if (j >= 2) {
								make_int(i, j - 2, taa, &a, CVEC(at_i->x), CVEC(at_j->x), &iout);
								double t3 = j * (j - 1) * t2;
								xin[idx + j + 60] = iout.x * t3;
								yin[idx + j + 60] = iout.y * t3;
								zin[idx + j + 60] = iout.z * t3;
							}
							else {
								xin[idx + j + 60] = 0.0;
								yin[idx + j + 60] = 0.0;
								zin[idx + j + 60] = 0.0;
							}
						}
					}
					for (size_t i = 0; i < count; i++) {
						size_t nx = shift_x[i];
						size_t ny = shift_y[i];
						size_t nz = shift_z[i];
						double xyz = xin[nx] * yin[ny] * zin[nz];
						double add = (xin[nx + 30] + xin[nx + 60]) * yin[ny] * zin[nz] +
							     (yin[ny + 30] + yin[ny + 60]) * xin[nx] * zin[nz] +
							     (zin[nz + 30] + zin[nz + 60]) * xin[nx] * yin[ny];
						sblk[i] = sblk[i] + dij[i] * xyz;
						tblk[i] = tblk[i] + dij[i] * (xyz * aj * ft[i] + add);
					}
				}
			}

			/* store integrals */
			for (size_t i = 0, idx = 0; i < count_i; i++) {
				size_t idx2 = (loc_i + i) * stride + loc_j;

				for (size_t j = 0; j < count_j; j++, idx++, idx2++) {
					s[idx2] = sblk[idx];
					t[idx2] = tblk[idx];
				}
			}
			loc_j += count_j;
		}}
		loc_i += count_i;
	}}
}

void efp_st_int_deriv(size_t n_atoms_i, const struct xr_atom *atoms_i,
		      size_t n_atoms_j, const struct xr_atom *atoms_j,
		      const vec_t *com_i, size_t size_i, size_t size_j,
		      six_t *ds, six_t *dt)
{
	static const size_t shift_x[] = { 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
					  3, 0, 0, 2, 2, 1, 0, 1, 0, 1 };

	static const size_t shift_y[] = { 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
					  0, 3, 0, 1, 0, 2, 2, 0, 1, 1 };

	static const size_t shift_z[] = { 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
					  0, 0, 3, 0, 1, 0, 1, 2, 2, 1 };

	double dij[100];
	double xs[5][6], ys[5][6], zs[5][6];
	double xt[5][4], yt[5][4], zt[5][4];
	double dxs[4][4], dys[4][4], dzs[4][4];
	double dxt[4][4], dyt[4][4], dzt[4][4];

	memset(ds, 0, size_i * size_j * sizeof(six_t));
	memset(dt, 0, size_i * size_j * sizeof(six_t));

	for (size_t iii = 0, loc_i = 0; iii < n_atoms_i; iii++) {
		const struct xr_atom *at_i = atoms_i + iii;

	/* shell i */
	for (size_t ii = 0; ii < at_i->n_shells; ii++) {
		const struct shell *sh_i = at_i->shells + ii;

		size_t type_i = get_shell_idx(sh_i->type);
		size_t start_i = get_shell_start(type_i);
		size_t end_i = get_shell_end(type_i);
		size_t sl_i = get_shell_sl(type_i);
		size_t count_i = end_i - start_i;

		for (size_t jjj = 0, loc_j = 0; jjj < n_atoms_j; jjj++) {
			const struct xr_atom *at_j = atoms_j + jjj;

		/* shell j */
		for (size_t jj = 0; jj < at_j->n_shells; jj++) {
			const struct shell *sh_j = at_j->shells + jj;

			size_t type_j = get_shell_idx(sh_j->type);
			size_t start_j = get_shell_start(type_j);
			size_t end_j = get_shell_end(type_j);
			size_t sl_j = get_shell_sl(type_j);
			size_t count_j = end_j - start_j;

			double rr = vec_dist_2(CVEC(at_i->x), CVEC(at_j->x));

			const double *coef_i = sh_i->coef;

			/* primitive i */
			for (size_t ig = 0; ig < sh_i->n_funcs; ig++) {
				double ai = *coef_i++;

				double con_i[20];
				set_coef(con_i, sh_i->type, coef_i);

				coef_i++;
				if (sh_i->type == 'L')
					coef_i++;

				const double *coef_j = sh_j->coef;

				/* primitive j */
				for (size_t jg = 0; jg < sh_j->n_funcs; jg++) {
					double aj = *coef_j++;

					double aa = 1.0 / (ai + aj);
					double tmp = ai * aj * rr * aa;

					if (tmp > int_tol) {
						coef_j++;
						if (sh_j->type == 'L')
							coef_j++;

						continue;
					}

					double con_j[20];
					set_coef(con_j, sh_j->type, coef_j);

					coef_j++;
					if (sh_j->type == 'L')
						coef_j++;

					double fac = exp(-tmp);

					for (size_t i = start_i, idx = 0; i < end_i; i++)
						for (size_t j = start_j; j < end_j; j++, idx++)
							dij[idx] = fac * con_i[i] * int_norm[i] * con_j[j] * int_norm[j];

					double taa = sqrt(aa);

					vec_t a = {
						(ai * at_i->x + aj * at_j->x) * aa,
						(ai * at_i->y + aj * at_j->y) * aa,
						(ai * at_i->z + aj * at_j->z) * aa
					};

					for (size_t i = 0; i < sl_i + 1; i++) {
						for (size_t j = 0; j < sl_j + 2; j++) {
							vec_t iout;
							make_int(i, j, taa, &a, CVEC(at_i->x), CVEC(at_j->x), &iout);
							xs[i][j] = iout.x * taa;
							ys[i][j] = iout.y * taa;
							zs[i][j] = iout.z * taa;
						}
					}

					double ai2 = 2.0 * ai;
					double aj2 = 2.0 * aj;

					for (size_t i = 0; i < sl_i + 1; i++) {
						xt[i][0] = (xs[i][0] - xs[i][2] * aj2) * aj;
						yt[i][0] = (ys[i][0] - ys[i][2] * aj2) * aj;
						zt[i][0] = (zs[i][0] - zs[i][2] * aj2) * aj;
					}

					if (sl_j > 1) {
						for (size_t i = 0; i < sl_i + 1; i++) {
							xt[i][1] = (xs[i][1] * 3.0 - xs[i][3] * aj2) * aj;
							yt[i][1] = (ys[i][1] * 3.0 - ys[i][3] * aj2) * aj;
							zt[i][1] = (zs[i][1] * 3.0 - zs[i][3] * aj2) * aj;
						}

						for (size_t j = 2; j < sl_j; j++) {
							for (size_t i = 0; i < sl_i + 1; i++) {
								size_t n1 = 2 * j + 1;
								size_t n2 = j * (j - 1) / 2;
								xt[i][j] = (xs[i][j] * n1 - xs[i][j + 2] * aj2) * aj - xs[i][j - 2] * n2;
								yt[i][j] = (ys[i][j] * n1 - ys[i][j + 2] * aj2) * aj - ys[i][j - 2] * n2;
								zt[i][j] = (zs[i][j] * n1 - zs[i][j + 2] * aj2) * aj - zs[i][j - 2] * n2;
							}
						}
					}

					for (size_t j = 0; j < sl_j; j++) {
						dxs[0][j] = xs[1][j] * ai2;
						dys[0][j] = ys[1][j] * ai2;
						dzs[0][j] = zs[1][j] * ai2;

						dxt[0][j] = xt[1][j] * ai2;
						dyt[0][j] = yt[1][j] * ai2;
						dzt[0][j] = zt[1][j] * ai2;
					}

					for (size_t i = 1; i < sl_i; i++) {
						for (size_t j = 0; j < sl_j; j++) {
							dxs[i][j] = xs[i + 1][j] * ai2 - xs[i - 1][j] * i;
							dys[i][j] = ys[i + 1][j] * ai2 - ys[i - 1][j] * i;
							dzs[i][j] = zs[i + 1][j] * ai2 - zs[i - 1][j] * i;

							dxt[i][j] = xt[i + 1][j] * ai2 - xt[i - 1][j] * i;
							dyt[i][j] = yt[i + 1][j] * ai2 - yt[i - 1][j] * i;
							dzt[i][j] = zt[i + 1][j] * ai2 - zt[i - 1][j] * i;
						}
					}

					for (size_t i = start_i, idx = 0; i < end_i; i++) {
						size_t ix = shift_x[i];
						size_t iy = shift_y[i];
						size_t iz = shift_z[i];

						for (size_t j = start_j; j < end_j; j++, idx++) {
							size_t jx = shift_x[j];
							size_t jy = shift_y[j];
							size_t jz = shift_z[j];

							double txs = dxs[ix][jx] * ys[iy][jy] * zs[iz][jz];
							double tys = xs[ix][jx] * dys[iy][jy] * zs[iz][jz];
							double tzs = xs[ix][jx] * ys[iy][jy] * dzs[iz][jz];

							double txt = dxt[ix][jx] * ys[iy][jy] * zs[iz][jz] +
								     dxs[ix][jx] * yt[iy][jy] * zs[iz][jz] +
								     dxs[ix][jx] * ys[iy][jy] * zt[iz][jz];
							double tyt = xt[ix][jx] * dys[iy][jy] * zs[iz][jz] +
								     xs[ix][jx] * dyt[iy][jy] * zs[iz][jz] +
								     xs[ix][jx] * dys[iy][jy] * zt[iz][jz];
							double tzt = xt[ix][jx] * ys[iy][jy] * dzs[iz][jz] +
								     xs[ix][jx] * yt[iy][jy] * dzs[iz][jz] +
								     xs[ix][jx] * ys[iy][jy] * dzt[iz][jz];

							size_t idx2 = (loc_i + i - start_i) * size_j + (loc_j + j - start_j);

							ds[idx2].x += txs * dij[idx];
							ds[idx2].y += tys * dij[idx];
							ds[idx2].z += tzs * dij[idx];
							ds[idx2].a += (tys * (at_i->z - com_i->z) - tzs * (at_i->y - com_i->y)) * dij[idx];
							ds[idx2].b += (tzs * (at_i->x - com_i->x) - txs * (at_i->z - com_i->z)) * dij[idx];
							ds[idx2].c += (txs * (at_i->y - com_i->y) - tys * (at_i->x - com_i->x)) * dij[idx];

							dt[idx2].x += txt * dij[idx];
							dt[idx2].y += tyt * dij[idx];
							dt[idx2].z += tzt * dij[idx];
							dt[idx2].a += (tyt * (at_i->z - com_i->z) - tzt * (at_i->y - com_i->y)) * dij[idx];
							dt[idx2].b += (tzt * (at_i->x - com_i->x) - txt * (at_i->z - com_i->z)) * dij[idx];
							dt[idx2].c += (txt * (at_i->y - com_i->y) - tyt * (at_i->x - com_i->x)) * dij[idx];
						}
					}
				}
			}
			loc_j += count_j;
		}}
		loc_i += count_i;
	}}
}
