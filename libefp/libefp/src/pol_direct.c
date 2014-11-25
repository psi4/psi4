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

#include <stdlib.h>

#include "clapack.h"
#include "compat.h"
#include "private.h"

double efp_get_pol_damp_tt(double, double, double);
enum efp_result efp_compute_id_direct(struct efp *);

static void
copy_matrix(double *dst, size_t n, size_t off_i, size_t off_j, const mat_t *m)
{
	dst[n * (3 * off_i + 0) + 3 * off_j + 0] = m->xx;
	dst[n * (3 * off_i + 0) + 3 * off_j + 1] = m->xy;
	dst[n * (3 * off_i + 0) + 3 * off_j + 2] = m->xz;
	dst[n * (3 * off_i + 1) + 3 * off_j + 0] = m->yx;
	dst[n * (3 * off_i + 1) + 3 * off_j + 1] = m->yy;
	dst[n * (3 * off_i + 1) + 3 * off_j + 2] = m->yz;
	dst[n * (3 * off_i + 2) + 3 * off_j + 0] = m->zx;
	dst[n * (3 * off_i + 2) + 3 * off_j + 1] = m->zy;
	dst[n * (3 * off_i + 2) + 3 * off_j + 2] = m->zz;
}

static void
transpose_matrix(double *m, size_t n)
{
	for (size_t i = 0; i < n; i++)
		for (size_t j = i + 1; j < n; j++) {
			double t = m[n * i + j];
			m[n * i + j] = m[n * j + i];
			m[n * j + i] = t;
		}
}

static mat_t
get_int_mat(const struct efp *efp, size_t i, size_t j, size_t ii, size_t jj)
{
	mat_t m;
	const struct frag *fr_i = efp->frags + i;
	const struct frag *fr_j = efp->frags + j;
	const struct polarizable_pt *pt_i = fr_i->polarizable_pts + ii;
	const struct polarizable_pt *pt_j = fr_j->polarizable_pts + jj;
	struct swf swf = efp_make_swf(efp, fr_i, fr_j);

	vec_t dr = {
		pt_j->x - pt_i->x - swf.cell.x,
		pt_j->y - pt_i->y - swf.cell.y,
		pt_j->z - pt_i->z - swf.cell.z
	};

	double p1 = 1.0;
	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;

	if (efp->opts.pol_damp == EFP_POL_DAMP_TT)
		p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp, fr_j->pol_damp);

	m.xx = swf.swf * p1 * (3.0 * dr.x * dr.x / r5 - 1.0 / r3);
	m.xy = swf.swf * p1 *  3.0 * dr.x * dr.y / r5;
	m.xz = swf.swf * p1 *  3.0 * dr.x * dr.z / r5;
	m.yx = swf.swf * p1 *  3.0 * dr.y * dr.x / r5;
	m.yy = swf.swf * p1 * (3.0 * dr.y * dr.y / r5 - 1.0 / r3);
	m.yz = swf.swf * p1 *  3.0 * dr.y * dr.z / r5;
	m.zx = swf.swf * p1 *  3.0 * dr.z * dr.x / r5;
	m.zy = swf.swf * p1 *  3.0 * dr.z * dr.y / r5;
	m.zz = swf.swf * p1 * (3.0 * dr.z * dr.z / r5 - 1.0 / r3);

	return (m);
}

static void
compute_lhs(const struct efp *efp, double *c, bool conj)
{
	size_t n = 3 * efp->n_polarizable_pts;

	for (size_t i = 0, offset_i = 0; i < efp->n_frag; i++) {
	for (size_t ii = 0; ii < efp->frags[i].n_polarizable_pts; ii++, offset_i++) {
	for (size_t j = 0, offset_j = 0; j < efp->n_frag; j++) {
	for (size_t jj = 0; jj < efp->frags[j].n_polarizable_pts; jj++, offset_j++) {
		if (i == j) {
			if (ii == jj)
				copy_matrix(c, n, offset_i, offset_j, &mat_identity);
			else
				copy_matrix(c, n, offset_i, offset_j, &mat_zero);

			continue;
		}

		const struct polarizable_pt *pt_i = efp->frags[i].polarizable_pts + ii;
		mat_t m = get_int_mat(efp, i, j, ii, jj);

		if (conj)
			m = mat_trans_mat(&pt_i->tensor, &m);
		else
			m = mat_mat(&pt_i->tensor, &m);

		mat_negate(&m);
		copy_matrix(c, n, offset_i, offset_j, &m);
	}}}}
}

static void
compute_rhs(const struct efp *efp, vec_t *id, bool conj)
{
	for (size_t i = 0, idx = 0; i < efp->n_frag; i++) {
		const struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++, idx++) {
			const struct polarizable_pt *pt = frag->polarizable_pts + j;
			vec_t field = vec_add(&pt->elec_field, &pt->elec_field_wf);

			if (conj)
				id[idx] = mat_trans_vec(&pt->tensor, &field);
			else
				id[idx] = mat_vec(&pt->tensor, &field);
		}
	}
}

enum efp_result
efp_compute_id_direct(struct efp *efp)
{
	double *c;
	size_t n;
	int *ipiv;
	enum efp_result res;

	n = 3 * efp->n_polarizable_pts;
	c = (double *)calloc(n * n, sizeof(double));
	ipiv = (int *)calloc(n, sizeof(int));

	if (c == NULL || ipiv == NULL) {
		res = EFP_RESULT_NO_MEMORY;
		goto error;
	}

	/* induced dipoles */
	compute_lhs(efp, c, false);
	compute_rhs(efp, efp->indip, false);
	transpose_matrix(c, n);

	if (efp_dgesv((int)n, 1, c, (int)n, ipiv, (double *)efp->indip, (int)n) != 0) {
		efp_log("dgesv: error solving for induced dipoles");
		res = EFP_RESULT_FATAL;
		goto error;
	}

	/* conjugate induced dipoles */
	compute_lhs(efp, c, true);
	compute_rhs(efp, efp->indipconj, true);
	transpose_matrix(c, n);

	if (efp_dgesv((int)n, 1, c, (int)n, ipiv, (double *)efp->indipconj, (int)n) != 0) {
		efp_log("dgesv: error solving for conjugate induced dipoles");
		res = EFP_RESULT_FATAL;
		goto error;
	}

	res = EFP_RESULT_SUCCESS;
error:
	free(c);
	free(ipiv);

	return (res);
}
