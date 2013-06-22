/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
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

#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "private.h"

#define INTEGRAL_THRESHOLD 1.0e-7

static inline int
fock_idx(int i, int j)
{
	return i < j ?
		j * (j + 1) / 2 + i :
		i * (i + 1) / 2 + j;
}

static inline double
valence(double n)
{
	if (n > 2)
		n -= 2;
	if (n > 8)
		n -= 8;
	return n;
}

static double
charge_penetration_energy(double s_ij, double r_ij)
{
	if (fabs(s_ij) < INTEGRAL_THRESHOLD)
		return 0.0;

	double ln_s = log(fabs(s_ij));

	return -s_ij * s_ij / r_ij / sqrt(-2.0 * ln_s);
}

static void
charge_penetration_grad(struct efp *efp, struct frag *fr_i, struct frag *fr_j,
			int lmo_i_idx, int lmo_j_idx, double s_ij, const six_t ds_ij,
			const struct swf *swf)
{
	if (fabs(s_ij) < INTEGRAL_THRESHOLD)
		return;

	const vec_t *ct_i = fr_i->lmo_centroids + lmo_i_idx;
	const vec_t *ct_j = fr_j->lmo_centroids + lmo_j_idx;

	vec_t dr = {
		ct_j->x - ct_i->x - swf->cell.x,
		ct_j->y - ct_i->y - swf->cell.y,
		ct_j->z - ct_i->z - swf->cell.z
	};

	double r_ij = vec_len(&dr);
	double ln_s = log(fabs(s_ij));

	double t1 = s_ij * s_ij / (r_ij * r_ij * r_ij) / sqrt(-2.0 * ln_s);
	double t2 = -s_ij / r_ij / sqrt(2.0) * (2.0 / sqrt(-ln_s) +
				0.5 / sqrt(-ln_s * ln_s * ln_s));

	vec_t force, torque_i, torque_j;

	force.x = (t2 * ds_ij.x - t1 * dr.x) * swf->swf;
	force.y = (t2 * ds_ij.y - t1 * dr.y) * swf->swf;
	force.z = (t2 * ds_ij.z - t1 * dr.z) * swf->swf;

	torque_i.x = swf->swf * (-t2 * ds_ij.a + t1 * (dr.y * (ct_i->z - fr_i->z) -
					dr.z * (ct_i->y - fr_i->y)));
	torque_i.y = swf->swf * (-t2 * ds_ij.b + t1 * (dr.z * (ct_i->x - fr_i->x) -
					dr.x * (ct_i->z - fr_i->z)));
	torque_i.z = swf->swf * (-t2 * ds_ij.c + t1 * (dr.x * (ct_i->y - fr_i->y) -
					dr.y * (ct_i->x - fr_i->x)));

	torque_j.x = torque_i.x + force.y * (fr_j->z - fr_i->z - swf->cell.z) -
				  force.z * (fr_j->y - fr_i->y - swf->cell.y);
	torque_j.y = torque_i.y + force.z * (fr_j->x - fr_i->x - swf->cell.x) -
				  force.x * (fr_j->z - fr_i->z - swf->cell.z);
	torque_j.z = torque_i.z + force.x * (fr_j->y - fr_i->y - swf->cell.y) -
				  force.y * (fr_j->x - fr_i->x - swf->cell.x);

	vec_atomic_add(&fr_i->force, &force);
	vec_atomic_add(&fr_i->torque, &torque_i);

	vec_atomic_sub(&fr_j->force, &force);
	vec_atomic_sub(&fr_j->torque, &torque_j);

	efp_add_stress(&swf->dr, &force, &efp->stress);
}

static void
transform_integrals(int n_lmo_i, int n_lmo_j,
		    int wf_size_i, int wf_size_j,
		    const double *wf_i, const double *wf_j,
		    const double *s, double *lmo_s)
{
	double tmp[n_lmo_i * wf_size_j];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		n_lmo_i, wf_size_j, wf_size_i, 1.0, wf_i, wf_size_i,
		s, wf_size_j, 0.0, tmp, wf_size_j);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		n_lmo_i, n_lmo_j, wf_size_j, 1.0, tmp, wf_size_j,
		wf_j, wf_size_j, 0.0, lmo_s, n_lmo_j);
}

static void
transform_integral_derivatives(int n_lmo_i, int n_lmo_j, int wf_size_i, int wf_size_j,
	const double *wf_i, const double *wf_j, const six_t *ds, six_t *lmo_ds)
{
	six_t tmp[n_lmo_i * wf_size_j];
	const six_t *p_ds;
	six_t *p_tmp, *p_lmo_ds;

	p_tmp = tmp;

	for (int i = 0; i < n_lmo_i; i++) {
		for (int j = 0; j < wf_size_j; j++, p_tmp++) {
			six_t sum = six_zero;
			p_ds = ds + j;

			for (int k = 0; k < wf_size_i; k++, p_ds += wf_size_j) {
				double w = wf_i[wf_size_i * i + k];

				sum.x += p_ds->x * w;
				sum.y += p_ds->y * w;
				sum.z += p_ds->z * w;
				sum.a += p_ds->a * w;
				sum.b += p_ds->b * w;
				sum.c += p_ds->c * w;
			}

			*p_tmp = sum;
		}
	}

	p_lmo_ds = lmo_ds;

	for (int i = 0; i < n_lmo_i; i++) {
		for (int j = 0; j < n_lmo_j; j++, p_lmo_ds++) {
			six_t sum = six_zero;
			p_tmp = tmp + i * wf_size_j;

			for (int k = 0; k < wf_size_j; k++, p_tmp++) {
				double w = wf_j[wf_size_j * j + k];

				sum.x += p_tmp->x * w;
				sum.y += p_tmp->y * w;
				sum.z += p_tmp->z * w;
				sum.a += p_tmp->a * w;
				sum.b += p_tmp->b * w;
				sum.c += p_tmp->c * w;
			}

			*p_lmo_ds = sum;
		}
	}
}

static void
add_six_vec(int el, int size, const double *vec, six_t *six)
{
	double *ptr = (double *)six + el;

	for (int i = 0; i < size; i++, vec++, ptr += 6)
		*ptr += *vec;
}

/*
 * Reference:
 *
 * Hui Li, Mark Gordon
 *
 * Gradients of the exchange-repulsion energy in the general effective fragment
 * potential method
 *
 * Theor. Chem. Acc. 115, 385 (2006)
 */
static void
lmo_lmo_xr_grad(struct efp *efp, struct frag *fr_i, struct frag *fr_j,
		int i, int j, const double *lmo_s, const double *lmo_t,
		const six_t *lmo_ds, const six_t *lmo_dt, const struct swf *swf)
{
	int ij = i * fr_j->n_lmo + j;

	const vec_t *ct_i = fr_i->lmo_centroids + i;
	const vec_t *ct_j = fr_j->lmo_centroids + j;

	vec_t dr = {
		ct_j->x - ct_i->x - swf->cell.x,
		ct_j->y - ct_i->y - swf->cell.y,
		ct_j->z - ct_i->z - swf->cell.z
	};

	double s_ij = lmo_s[ij];
	double t_ij = lmo_t[ij];
	double r_ij = vec_len(&dr);
	double r_ij3 = r_ij * r_ij * r_ij;

	six_t ds_ij = lmo_ds[ij];
	six_t dt_ij = lmo_dt[ij];

	double t1, t2;
	vec_t force = vec_zero, torque_i = vec_zero;

	/* first part */
	if (fabs(s_ij) > INTEGRAL_THRESHOLD) {
		double ln_s = log(fabs(s_ij));

		t1 = s_ij / r_ij * (-sqrt(-2.0 / PI / ln_s) + 4.0 * sqrt(-2.0 / PI * ln_s));
		t2 = 2.0 * sqrt(-2.0 / PI * ln_s) * s_ij * s_ij / r_ij3;

		force.x += -t1 * ds_ij.x - t2 * dr.x;
		force.y += -t1 * ds_ij.y - t2 * dr.y;
		force.z += -t1 * ds_ij.z - t2 * dr.z;

		torque_i.x += t1 * ds_ij.a + t2 * (dr.y * (ct_i->z - fr_i->z) -
						   dr.z * (ct_i->y - fr_i->y));
		torque_i.y += t1 * ds_ij.b + t2 * (dr.z * (ct_i->x - fr_i->x) -
						   dr.x * (ct_i->z - fr_i->z));
		torque_i.z += t1 * ds_ij.c + t2 * (dr.x * (ct_i->y - fr_i->y) -
						   dr.y * (ct_i->x - fr_i->x));
	}

	/* second part */
	double fij = 0.0, fji = 0.0;
	six_t dfij = six_zero, dfji = six_zero;

	for (int k = 0; k < fr_i->n_lmo; k++) {
		double fe = fr_i->xr_fock_mat[fock_idx(i, k)];
		const six_t *ds = lmo_ds + k * fr_j->n_lmo + j;

		fij += fe * lmo_s[k * fr_j->n_lmo + j];

		dfij.x += fe * ds->x;
		dfij.y += fe * ds->y;
		dfij.z += fe * ds->z;
		dfij.a += fe * ds->a;
		dfij.b += fe * ds->b;
		dfij.c += fe * ds->c;
	}

	for (int l = 0; l < fr_j->n_lmo; l++) {
		double fe = fr_j->xr_fock_mat[fock_idx(j, l)];
		const six_t *ds = lmo_ds + i * fr_j->n_lmo + l;

		fji += fe * lmo_s[i * fr_j->n_lmo + l];

		dfji.x += fe * ds->x;
		dfji.y += fe * ds->y;
		dfji.z += fe * ds->z;
		dfji.a += fe * ds->a;
		dfji.b += fe * ds->b;
		dfji.c += fe * ds->c;
	}

	t1 = (fij + fji - 2.0 * t_ij);

	force.x += -t1 * ds_ij.x - s_ij * (dfij.x + dfji.x - 2.0 * dt_ij.x);
	force.y += -t1 * ds_ij.y - s_ij * (dfij.y + dfji.y - 2.0 * dt_ij.y);
	force.z += -t1 * ds_ij.z - s_ij * (dfij.z + dfji.z - 2.0 * dt_ij.z);

	torque_i.x += t1 * ds_ij.a + s_ij * (dfij.a + dfji.a - 2.0 * dt_ij.a);
	torque_i.y += t1 * ds_ij.b + s_ij * (dfij.b + dfji.b - 2.0 * dt_ij.b);
	torque_i.z += t1 * ds_ij.c + s_ij * (dfij.c + dfji.c - 2.0 * dt_ij.c);

	/* third part */
	double vib = 0.0, vja = 0.0;
	six_t dvib = six_zero, dvja = six_zero;

	for (int l = 0; l < fr_j->n_atoms; l++) {
		struct efp_atom *at_j = fr_j->atoms + l;

		vec_t dr_a = {
			at_j->x - ct_i->x - swf->cell.x,
			at_j->y - ct_i->y - swf->cell.y,
			at_j->z - ct_i->z - swf->cell.z
		};

		double r = vec_len(&dr_a);
		double tmp = valence(at_j->znuc) / (r * r * r);

		vib -= valence(at_j->znuc) / r;

		dvib.x -= tmp * dr_a.x;
		dvib.y -= tmp * dr_a.y;
		dvib.z -= tmp * dr_a.z;

		dvib.a -= tmp * (dr_a.y * (ct_i->z - fr_i->z) -
				 dr_a.z * (ct_i->y - fr_i->y));
		dvib.b -= tmp * (dr_a.z * (ct_i->x - fr_i->x) -
				 dr_a.x * (ct_i->z - fr_i->z));
		dvib.c -= tmp * (dr_a.x * (ct_i->y - fr_i->y) -
				 dr_a.y * (ct_i->x - fr_i->x));
	}

	for (int l = 0; l < fr_j->n_lmo; l++) {
		vec_t *ct_jj = fr_j->lmo_centroids + l;

		vec_t dr_a = {
			ct_jj->x - ct_i->x - swf->cell.x,
			ct_jj->y - ct_i->y - swf->cell.y,
			ct_jj->z - ct_i->z - swf->cell.z
		};

		double r = vec_len(&dr_a);
		double tmp = 2.0 / (r * r * r);

		vib += 2.0 / r;

		dvib.x += tmp * dr_a.x;
		dvib.y += tmp * dr_a.y;
		dvib.z += tmp * dr_a.z;

		dvib.a += tmp * (dr_a.y * (ct_i->z - fr_i->z) -
				 dr_a.z * (ct_i->y - fr_i->y));
		dvib.b += tmp * (dr_a.z * (ct_i->x - fr_i->x) -
				 dr_a.x * (ct_i->z - fr_i->z));
		dvib.c += tmp * (dr_a.x * (ct_i->y - fr_i->y) -
				 dr_a.y * (ct_i->x - fr_i->x));
	}

	for (int k = 0; k < fr_i->n_atoms; k++) {
		struct efp_atom *at_i = fr_i->atoms + k;

		vec_t dr_a = {
			ct_j->x - at_i->x - swf->cell.x,
			ct_j->y - at_i->y - swf->cell.y,
			ct_j->z - at_i->z - swf->cell.z
		};

		double r = vec_len(&dr_a);
		double tmp = valence(at_i->znuc) / (r * r * r);

		vja -= valence(at_i->znuc) / r;

		dvja.x -= tmp * dr_a.x;
		dvja.y -= tmp * dr_a.y;
		dvja.z -= tmp * dr_a.z;

		dvja.a -= tmp * (dr_a.y * (at_i->z - fr_i->z) -
				 dr_a.z * (at_i->y - fr_i->y));
		dvja.b -= tmp * (dr_a.z * (at_i->x - fr_i->x) -
				 dr_a.x * (at_i->z - fr_i->z));
		dvja.c -= tmp * (dr_a.x * (at_i->y - fr_i->y) -
				 dr_a.y * (at_i->x - fr_i->x));
	}

	for (int k = 0; k < fr_i->n_lmo; k++) {
		vec_t *ct_ii = fr_i->lmo_centroids + k;

		vec_t dr_a = {
			ct_j->x - ct_ii->x - swf->cell.x,
			ct_j->y - ct_ii->y - swf->cell.y,
			ct_j->z - ct_ii->z - swf->cell.z
		};

		double r = vec_len(&dr_a);
		double tmp = 2.0 / (r * r * r);

		vja += 2.0 / r;

		dvja.x += tmp * dr_a.x;
		dvja.y += tmp * dr_a.y;
		dvja.z += tmp * dr_a.z;

		dvja.a += tmp * (dr_a.y * (ct_ii->z - fr_i->z) -
				 dr_a.z * (ct_ii->y - fr_i->y));
		dvja.b += tmp * (dr_a.z * (ct_ii->x - fr_i->x) -
				 dr_a.x * (ct_ii->z - fr_i->z));
		dvja.c += tmp * (dr_a.x * (ct_ii->y - fr_i->y) -
				 dr_a.y * (ct_ii->x - fr_i->x));
	}

	t1 = 2.0 * s_ij * (vib + vja - 1.0 / r_ij);

	force.x += t1 * ds_ij.x + s_ij * s_ij * (dvib.x + dvja.x - dr.x / r_ij3);
	force.y += t1 * ds_ij.y + s_ij * s_ij * (dvib.y + dvja.y - dr.y / r_ij3);
	force.z += t1 * ds_ij.z + s_ij * s_ij * (dvib.z + dvja.z - dr.z / r_ij3);

	torque_i.x += -t1 * ds_ij.a - s_ij * s_ij * (dvib.a + dvja.a -
			(dr.y * (ct_i->z - fr_i->z) -
			 dr.z * (ct_i->y - fr_i->y)) / r_ij3);
	torque_i.y += -t1 * ds_ij.b - s_ij * s_ij * (dvib.b + dvja.b -
			(dr.z * (ct_i->x - fr_i->x) -
			 dr.x * (ct_i->z - fr_i->z)) / r_ij3);
	torque_i.z += -t1 * ds_ij.c - s_ij * s_ij * (dvib.c + dvja.c -
			(dr.x * (ct_i->y - fr_i->y) -
			 dr.y * (ct_i->x - fr_i->x)) / r_ij3);

	force.x *= 2.0 * swf->swf;
	force.y *= 2.0 * swf->swf;
	force.z *= 2.0 * swf->swf;

	torque_i.x *= 2.0 * swf->swf;
	torque_i.y *= 2.0 * swf->swf;
	torque_i.z *= 2.0 * swf->swf;

	vec_t torque_j = {
		torque_i.x + force.y * (fr_j->z - fr_i->z - swf->cell.z) -
			force.z * (fr_j->y - fr_i->y - swf->cell.y),
		torque_i.y + force.z * (fr_j->x - fr_i->x - swf->cell.x) -
			force.x * (fr_j->z - fr_i->z - swf->cell.z),
		torque_i.z + force.x * (fr_j->y - fr_i->y - swf->cell.y) -
			force.y * (fr_j->x - fr_i->x - swf->cell.x)
	};

	vec_atomic_add(&fr_i->force, &force);
	vec_atomic_add(&fr_i->torque, &torque_i);

	vec_atomic_sub(&fr_j->force, &force);
	vec_atomic_sub(&fr_j->torque, &torque_j);

	efp_add_stress(&swf->dr, &force, &efp->stress);
}

static double
lmo_lmo_xr_energy(struct frag *fr_i, struct frag *fr_j, int i, int j,
		  const double *lmo_s, const double *lmo_t, const struct swf *swf)
{
	double s_ij = lmo_s[i * fr_j->n_lmo + j];
	double t_ij = lmo_t[i * fr_j->n_lmo + j];

	const vec_t *ct_i = fr_i->lmo_centroids + i;
	const vec_t *ct_j = fr_j->lmo_centroids + j;

	vec_t dr = {
		ct_j->x - ct_i->x - swf->cell.x,
		ct_j->y - ct_i->y - swf->cell.y,
		ct_j->z - ct_i->z - swf->cell.z
	};

	double r_ij = vec_len(&dr);
	double exr = 0.0;

	/* xr - first part */
	if (fabs(s_ij) > INTEGRAL_THRESHOLD)
		exr += -2.0 * sqrt(-2.0 * log(fabs(s_ij)) / PI) * s_ij * s_ij / r_ij;

	/* xr - second part */
	for (int k = 0; k < fr_i->n_lmo; k++) {
		exr -= s_ij * lmo_s[k * fr_j->n_lmo + j] *
				fr_i->xr_fock_mat[fock_idx(i, k)];
	}
	for (int l = 0; l < fr_j->n_lmo; l++) {
		exr -= s_ij * lmo_s[i * fr_j->n_lmo + l] *
				fr_j->xr_fock_mat[fock_idx(j, l)];
	}
	exr += 2.0 * s_ij * t_ij;

	/* xr - third part */
	for (int jj = 0; jj < fr_j->n_atoms; jj++) {
		struct efp_atom *at_j = fr_j->atoms + jj;

		vec_t dr_a = {
			at_j->x - ct_i->x - swf->cell.x,
			at_j->y - ct_i->y - swf->cell.y,
			at_j->z - ct_i->z - swf->cell.z
		};

		double r = vec_len(&dr_a);

		exr -= s_ij * s_ij * valence(at_j->znuc) / r;
	}
	for (int jj = 0; jj < fr_j->n_lmo; jj++) {
		const vec_t *ct_jj = fr_j->lmo_centroids + jj;

		vec_t dr_a = {
			ct_jj->x - ct_i->x - swf->cell.x,
			ct_jj->y - ct_i->y - swf->cell.y,
			ct_jj->z - ct_i->z - swf->cell.z
		};

		double r = vec_len(&dr_a);

		exr += 2.0 * s_ij * s_ij / r;
	}
	for (int ii = 0; ii < fr_i->n_atoms; ii++) {
		struct efp_atom *at_i = fr_i->atoms + ii;

		vec_t dr_a = {
			ct_j->x - at_i->x - swf->cell.x,
			ct_j->y - at_i->y - swf->cell.y,
			ct_j->z - at_i->z - swf->cell.z
		};

		double r = vec_len(&dr_a);

		exr -= s_ij * s_ij * valence(at_i->znuc) / r;
	}
	for (int ii = 0; ii < fr_i->n_lmo; ii++) {
		const vec_t *ct_ii = fr_i->lmo_centroids + ii;

		vec_t dr_a = {
			ct_j->x - ct_ii->x - swf->cell.x,
			ct_j->y - ct_ii->y - swf->cell.y,
			ct_j->z - ct_ii->z - swf->cell.z
		};

		double r = vec_len(&dr_a);

		exr += 2.0 * s_ij * s_ij / r;
	}
	exr -= s_ij * s_ij / r_ij;

	return 2.0 * exr;
}

static void
frag_frag_xr(struct efp *efp, int frag_i, int frag_j, int overlap_idx,
	     double *exr_out, double *ecp_out)
{
	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	double *s = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(double));
	double *t = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(double));

	double lmo_s[fr_i->n_lmo * fr_j->n_lmo];
	double lmo_t[fr_i->n_lmo * fr_j->n_lmo];

	struct swf swf = efp_make_swf(efp, fr_i, fr_j);
	struct shell shells_j[fr_j->n_xr_shells];

	for (int j = 0; j < fr_j->n_xr_shells; j++) {
		shells_j[j] = fr_j->xr_shells[j];

		shells_j[j].x -= swf.cell.x;
		shells_j[j].y -= swf.cell.y;
		shells_j[j].z -= swf.cell.z;
	}

	efp_st_int(fr_i->n_xr_shells, fr_i->xr_shells,
		   fr_j->n_xr_shells, shells_j,
		   fr_j->xr_wf_size, s, t);

	transform_integrals(fr_i->n_lmo, fr_j->n_lmo,
			    fr_i->xr_wf_size, fr_j->xr_wf_size,
			    fr_i->xr_wf, fr_j->xr_wf,
			    s, lmo_s);
	transform_integrals(fr_i->n_lmo, fr_j->n_lmo,
			    fr_i->xr_wf_size, fr_j->xr_wf_size,
			    fr_i->xr_wf, fr_j->xr_wf,
			    t, lmo_t);

	double exr = 0.0;
	double ecp = 0.0;

	for (int i = 0, idx = 0; i < fr_i->n_lmo; i++) {
		for (int j = 0; j < fr_j->n_lmo; j++, idx++) {
			double s_ij = lmo_s[i * fr_j->n_lmo + j];

			vec_t dr = {
				fr_j->lmo_centroids[j].x -
					fr_i->lmo_centroids[i].x - swf.cell.x,
				fr_j->lmo_centroids[j].y -
					fr_i->lmo_centroids[i].y - swf.cell.y,
				fr_j->lmo_centroids[j].z -
					fr_i->lmo_centroids[i].z - swf.cell.z
			};

			double r_ij = vec_len(&dr);

			if ((efp->opts.terms & EFP_TERM_DISP) &&
			    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP))
				fr_i->overlap_int[overlap_idx + idx] = s_ij;

			if ((efp->opts.terms & EFP_TERM_ELEC) &&
			    (efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP))
				ecp += charge_penetration_energy(s_ij, r_ij);

			if (efp->opts.terms & EFP_TERM_XR)
				exr += lmo_lmo_xr_energy(fr_i, fr_j, i, j,
							 lmo_s, lmo_t, &swf);
		}
	}

	*exr_out = exr * swf.swf;
	*ecp_out = ecp * swf.swf;

	if (!efp->do_gradient) {
		free(s);
		free(t);
		return;
	}

	/* compute gradient */

	six_t *ds = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(six_t));
	six_t *dt = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(six_t));

	six_t *lmo_ds = malloc(fr_i->n_lmo * fr_j->n_lmo * sizeof(six_t));
	six_t *lmo_dt = malloc(fr_i->n_lmo * fr_j->n_lmo * sizeof(six_t));

	efp_st_int_deriv(fr_i->n_xr_shells, fr_i->xr_shells,
			 fr_j->n_xr_shells, shells_j,
			 VEC(fr_i->x), fr_i->xr_wf_size, fr_j->xr_wf_size,
			 ds, dt);

	transform_integral_derivatives(fr_i->n_lmo, fr_j->n_lmo,
				       fr_i->xr_wf_size, fr_j->xr_wf_size,
				       fr_i->xr_wf, fr_j->xr_wf,
				       ds, lmo_ds);
	transform_integral_derivatives(fr_i->n_lmo, fr_j->n_lmo,
				       fr_i->xr_wf_size, fr_j->xr_wf_size,
				       fr_i->xr_wf, fr_j->xr_wf,
				       dt, lmo_dt);

	double lmo_tmp[fr_i->n_lmo * fr_j->n_lmo];

	for (int a = 0; a < 3; a++) {
		transform_integrals(fr_i->n_lmo, fr_j->n_lmo,
				    fr_i->xr_wf_size, fr_j->xr_wf_size,
				    fr_i->xr_wf_deriv[a], fr_j->xr_wf,
				    s, lmo_tmp);
		add_six_vec(3 + a, fr_i->n_lmo * fr_j->n_lmo, lmo_tmp, lmo_ds);

		transform_integrals(fr_i->n_lmo, fr_j->n_lmo,
				    fr_i->xr_wf_size, fr_j->xr_wf_size,
				    fr_i->xr_wf_deriv[a], fr_j->xr_wf,
				    t, lmo_tmp);
		add_six_vec(3 + a, fr_i->n_lmo * fr_j->n_lmo, lmo_tmp, lmo_dt);
	}

	for (int i = 0, idx = 0; i < fr_i->n_lmo; i++) {
		for (int j = 0; j < fr_j->n_lmo; j++, idx++) {
			int ij = i * fr_j->n_lmo + j;

			if ((efp->opts.terms & EFP_TERM_DISP) &&
			    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP))
				fr_i->overlap_int_deriv[overlap_idx + idx] = lmo_ds[ij];

			if ((efp->opts.terms & EFP_TERM_ELEC) &&
			    (efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP))
				charge_penetration_grad(efp, fr_i, fr_j, i, j,
							lmo_s[ij], lmo_ds[ij], &swf);

			if (efp->opts.terms & EFP_TERM_XR)
				lmo_lmo_xr_grad(efp, fr_i, fr_j, i, j, lmo_s, lmo_t,
						lmo_ds, lmo_dt, &swf);
		}
	}

	vec_t force = {
		swf.dswf.x * (exr + ecp),
		swf.dswf.y * (exr + ecp),
		swf.dswf.z * (exr + ecp)
	};

	vec_atomic_add(&fr_i->force, &force);
	vec_atomic_sub(&fr_j->force, &force);
	efp_add_stress(&swf.dr, &force, &efp->stress);

	free(s), free(ds), free(lmo_ds);
	free(t), free(dt), free(lmo_dt);
}

enum efp_result
efp_compute_xr(struct efp *efp)
{
	int do_xr = (efp->opts.terms & EFP_TERM_XR);
	int do_cp = (efp->opts.terms & EFP_TERM_ELEC) &&
		    (efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP);
	int do_dd = (efp->opts.terms & EFP_TERM_DISP) &&
		    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP);

	if (!do_xr && !do_cp && !do_dd)
		return EFP_RESULT_SUCCESS;

	double exr = 0.0;
	double ecp = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(+:exr,ecp)
#endif
	for (int i = 0; i < efp->n_frag; i++) {
		for (int j = i + 1, idx = 0; j < efp->n_frag; j++) {
			int n_lmo_ij = efp->frags[i].n_lmo * efp->frags[j].n_lmo;

			if (!efp_skip_frag_pair(efp, i, j)) {
				double exr_out, ecp_out;

				frag_frag_xr(efp, i, j, idx, &exr_out, &ecp_out);

				exr += exr_out;
				ecp += ecp_out;
			}
			else {
				/* also skips integrals */
				if ((efp->opts.terms & EFP_TERM_DISP) &&
				    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP)) {
					memset(efp->frags[i].overlap_int + idx, 0,
								n_lmo_ij * sizeof(double));
					memset(efp->frags[i].overlap_int_deriv + idx, 0,
								n_lmo_ij * sizeof(six_t));
				}
			}

			idx += n_lmo_ij;
		}
	}

	efp->energy.exchange_repulsion = exr;
	efp->energy.charge_penetration = ecp;

	return EFP_RESULT_SUCCESS;
}

static inline int
func_d_idx(int a, int b)
{
	/* order in which GAMESS stores D functions */
	enum { xx = 0, yy, zz, xy, xz, yz };

	static const int idx[] = {
		xx, xy, xz, xy, yy, yz, xz, yz, zz
	};

	return idx[3 * a + b];
}

static inline int
func_f_idx(int a, int b, int c)
{
	/* order in which GAMESS stores F functions */
	enum { xxx = 0, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz };

	static const int idx[] = {
		xxx, xxy, xxz, xxy, xyy, xyz, xxz, xyz, xzz,
		xxy, xyy, xyz, xyy, yyy, yyz, xyz, yyz, yzz,
		xxz, xyz, xzz, xyz, yyz, yzz, xzz, yzz, zzz
	};

	return idx[9 * a + 3 * b + c];
}

static void
rotate_func_d(const mat_t *rotmat, const double *in, double *out)
{
	const double norm = sqrt(3.0) / 2.0;

	double full_in[9], full_out[9];

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			full_in[3 * a + b] = in[func_d_idx(a, b)];
			if (a != b)
				full_in[3 * a + b] *= norm;
		}
	}

	efp_rotate_t2(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			if (a != b)
				full_out[3 * a + b] /= norm;
			out[func_d_idx(a, b)] = full_out[3 * a + b];
		}
	}
}

static void
rotate_func_f(const mat_t *rotmat, const double *in, double *out)
{
	const double norm1 = sqrt(5.0) / 3.0;
	const double norm2 = sqrt(3.0) / 2.0;

	double full_in[27], full_out[27];

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int full_idx = 9 * a + 3 * b + c;
				int in_idx = func_f_idx(a, b, c);

				full_in[full_idx] = in[in_idx];

				if (a != b || a != c)
					full_in[full_idx] *= norm1;

				if (a != b && a != c && b != c)
					full_in[full_idx] *= norm2;
			}

	efp_rotate_t3(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int full_idx = 9 * a + 3 * b + c;
				int out_idx = func_f_idx(a, b, c);

				if (a != b || a != c)
					full_out[full_idx] /= norm1;

				if (a != b && a != c && b != c)
					full_out[full_idx] /= norm2;

				out[out_idx] = full_out[full_idx];
			}
}

static void
coef_deriv_p(int axis, const double *coef, double *der)
{
	switch (axis) {
	case 0:
		der[0] = 0.0;
		der[1] = coef[2];
		der[2] = -coef[1];
		break;
	case 1:
		der[0] = -coef[2];
		der[1] = 0.0;
		der[2] = coef[0];
		break;
	case 2:
		der[0] = coef[1];
		der[1] = -coef[0];
		der[2] = 0.0;
		break;
	};
}

static void
coef_deriv_d(int axis, const double *coef, double *der)
{
	const double sqrt3 = sqrt(3.0);

	switch (axis) {
	case 0:
		der[0] = 0.0;
		der[1] = sqrt3 * coef[5];
		der[2] = -sqrt3 * coef[5];
		der[3] = coef[4];
		der[4] = -coef[3];
		der[5] = 2.0 / sqrt3 * (coef[2] - coef[1]);
		break;
	case 1:
		der[0] = -sqrt3 * coef[4];
		der[1] = 0.0;
		der[2] = sqrt3 * coef[4];
		der[3] = -coef[5];
		der[4] = 2.0 / sqrt3 * (coef[0] - coef[2]);
		der[5] = coef[3];
		break;
	case 2:
		der[0] = sqrt3 * coef[3];
		der[1] = -sqrt3 * coef[3];
		der[2] = 0.0;
		der[3] = 2.0 / sqrt3 * (coef[1] - coef[0]);
		der[4] = coef[5];
		der[5] = -coef[4];
		break;
	};
}

static void
coef_deriv_f(int axis, const double *coef, double *der)
{
	const double sqrt3 = sqrt(3.0);
	const double sqrt5 = sqrt(5.0);

	switch (axis) {
	case 0:
		der[0] = 0.0;
		der[1] = sqrt5 * coef[6];
		der[2] = -sqrt5 * coef[8];
		der[3] = coef[4];
		der[4] = -coef[3];
		der[5] = sqrt3 * coef[9];
		der[6] = -3.0 / sqrt5 * coef[1] + 2.0 * coef[8];
		der[7] = -sqrt3 * coef[9];
		der[8] = 3.0 / sqrt5 * coef[2] - 2.0 * coef[6];
		der[9] = 2.0 / sqrt3 * (coef[7] - coef[5]);
		break;
	case 1:
		der[0] = -sqrt5 * coef[4];
		der[1] = 0.0;
		der[2] = sqrt5 * coef[7];
		der[3] = -sqrt3 * coef[9];
		der[4] = 3.0 / sqrt5 * coef[0] - 2.0 * coef[7];
		der[5] = -coef[6];
		der[6] = coef[5];
		der[7] = -3.0 / sqrt5 * coef[2] + 2.0 * coef[4];
		der[8] = sqrt3 * coef[9];
		der[9] = 2.0 / sqrt3 * (coef[3] - coef[8]);
		break;
	case 2:
		der[0] = sqrt5 * coef[3];
		der[1] = -sqrt5 * coef[5];
		der[2] = 0.0;
		der[3] = -3.0 / sqrt5 * coef[0] + 2.0 * coef[5];
		der[4] = sqrt3 * coef[9];
		der[5] = 3.0 / sqrt5 * coef[1] - 2.0 * coef[3];
		der[6] = -sqrt3 * coef[9];
		der[7] = coef[8];
		der[8] = -coef[7];
		der[9] = 2.0 / sqrt3 * (coef[6] - coef[4]);
		break;
	};
}

void
efp_update_xr(struct frag *frag)
{
	const mat_t *rotmat = &frag->rotmat;

	/* update LMO centroids */
	for (int i = 0; i < frag->n_lmo; i++)
		efp_move_pt(CVEC(frag->x), rotmat, frag->lib->lmo_centroids + i,
				frag->lmo_centroids + i);

	/* update shells */
	for (int i = 0; i < frag->n_xr_shells; i++)
		efp_move_pt(CVEC(frag->x), rotmat, CVEC(frag->lib->xr_shells[i].x),
				VEC(frag->xr_shells[i].x));

	/* rotate wavefunction */
	for (int k = 0; k < frag->n_lmo; k++) {
		double *deriv[3];

		for (int a = 0; a < 3; a++)
			deriv[a] = frag->xr_wf_deriv[a] + k * frag->xr_wf_size;

		const double *in = frag->lib->xr_wf + k * frag->xr_wf_size;
		double *out = frag->xr_wf + k * frag->xr_wf_size;

		for (int i = 0, func = 0; i < frag->n_xr_shells; i++) {
			switch (frag->xr_shells[i].type) {
			case 'S':
				func++;
				break;
			case 'L':
				func++;
				/* fall through */
			case 'P': {
				vec_t r = mat_vec(rotmat, (const vec_t *)(in + func));

				out[func + 0] = r.x;
				out[func + 1] = r.y;
				out[func + 2] = r.z;

				for (int a = 0; a < 3; a++)
					coef_deriv_p(a, out + func, deriv[a] + func);

				func += 3;
				break;
			}
			case 'D':
				rotate_func_d(rotmat, in + func, out + func);

				for (int a = 0; a < 3; a++)
					coef_deriv_d(a, out + func, deriv[a] + func);

				func += 6;
				break;
			case 'F':
				rotate_func_f(rotmat, in + func, out + func);

				for (int a = 0; a < 3; a++)
					coef_deriv_f(a, out + func, deriv[a] + func);

				func += 10;
				break;
			}
		}
	}
}
