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

#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "efp_private.h"
#include "disp.h"

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

static inline double
get_disp_damp(double s_ij)
{
	if (fabs(s_ij) < 1.0e-6)
		return 0.0;

	double ln_s = log(fabs(s_ij));

	return 1.0 - s_ij * s_ij * (1.0 - 2.0 * ln_s + 2.0 * ln_s * ln_s);
}

static inline void
calc_disp_damp(struct efp *efp, int frag_i, int frag_j,
	       int i, int j, double s_ij)
{
	double damp = get_disp_damp(s_ij);

	efp->disp_damp_overlap[
		disp_damp_overlap_idx(efp, frag_i, frag_j, i, j)] = damp;
	efp->disp_damp_overlap[
		disp_damp_overlap_idx(efp, frag_j, frag_i, j, i)] = damp;
}

static inline double
get_charge_penetration(double s_ij, double r_ij)
{
	double ln_s = log(fabs(s_ij));
	return -2.0 * s_ij * s_ij / r_ij / sqrt(-2.0 * ln_s);
}

static void
frag_frag_xr_grad(struct efp *efp, int frag_i, int frag_j)
{
	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	int size = fr_i->xr_wf_size * fr_j->xr_wf_size;

	double *ds = malloc(6 * size * sizeof(double));
	double *dt = malloc(6 * size * sizeof(double));

	efp_st_int_deriv(fr_i->n_xr_shells, fr_i->xr_shells,
			 fr_j->n_xr_shells, fr_j->xr_shells,
			 VEC(fr_i->x), fr_i->xr_wf_size, fr_j->xr_wf_size,
			 ds, dt);

//	double *ds_x = ds + 0 * size;
//	double *ds_y = ds + 1 * size;
//	double *ds_z = ds + 2 * size;
//	double *ds_a = ds + 3 * size;
//	double *ds_b = ds + 4 * size;
//	double *ds_c = ds + 5 * size;
//
//	double *dt_x = dt + 0 * size;
//	double *dt_y = dt + 1 * size;
//	double *dt_z = dt + 2 * size;
//	double *dt_a = dt + 3 * size;
//	double *dt_b = dt + 4 * size;
//	double *dt_c = dt + 5 * size;

	free(ds);
	free(dt);
}

static void
frag_frag_xr(struct efp *efp, int frag_i, int frag_j,
	     double *exr_out, double *ecp_out)
{
	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	double s[fr_i->xr_wf_size * fr_j->xr_wf_size];
	double t[fr_i->xr_wf_size * fr_j->xr_wf_size];

	double lmo_s[fr_i->n_lmo][fr_j->n_lmo];
	double lmo_t[fr_i->n_lmo][fr_j->n_lmo];

	double tmp[fr_i->n_lmo * fr_j->xr_wf_size];

	/* compute S and T integrals */
	efp_st_int(fr_i->n_xr_shells, fr_i->xr_shells,
		   fr_j->n_xr_shells, fr_j->xr_shells,
		   fr_j->xr_wf_size, s, t);

	/* lmo_s = wf_i * s * wf_j(t) */
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		fr_i->n_lmo, fr_j->xr_wf_size, fr_i->xr_wf_size, 1.0,
		fr_i->xr_wf, fr_i->xr_wf_size, s, fr_j->xr_wf_size, 0.0,
		tmp, fr_j->xr_wf_size);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		fr_i->n_lmo, fr_j->n_lmo, fr_j->xr_wf_size, 1.0,
		tmp, fr_j->xr_wf_size, fr_j->xr_wf, fr_j->xr_wf_size, 0.0,
		&lmo_s[0][0], fr_j->n_lmo);

	/* lmo_t = wf_i * t * wf_j(t) */
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		fr_i->n_lmo, fr_j->xr_wf_size, fr_i->xr_wf_size, 1.0,
		fr_i->xr_wf, fr_i->xr_wf_size, t, fr_j->xr_wf_size, 0.0,
		tmp, fr_j->xr_wf_size);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		fr_i->n_lmo, fr_j->n_lmo, fr_j->xr_wf_size, 1.0,
		tmp, fr_j->xr_wf_size, fr_j->xr_wf, fr_j->xr_wf_size, 0.0,
		&lmo_t[0][0], fr_j->n_lmo);

	double exr = 0.0;
	double ecp = 0.0;

	for (int i = 0; i < fr_i->n_lmo; i++) {
		for (int j = 0; j < fr_j->n_lmo; j++) {
			double s_ij = lmo_s[i][j];
			double t_ij = lmo_t[i][j];
			double r_ij = vec_dist(fr_i->lmo_centroids + i,
					       fr_j->lmo_centroids + j);

			if ((efp->opts.terms & EFP_TERM_DISP) &&
			    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP)) {
				calc_disp_damp(efp, frag_i, frag_j, i, j, s_ij);
			}

			if ((efp->opts.terms & EFP_TERM_ELEC) &&
			    (efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP)) {
				ecp += get_charge_penetration(s_ij, r_ij);
			}

			/* xr - first part */
			if (fabs(s_ij) > 1.0e-6)
				exr += -2.0 * sqrt(-2.0 * log(fabs(s_ij)) / PI)
							* s_ij * s_ij / r_ij;

			/* xr - second part */
			for (int k = 0; k < fr_i->n_lmo; k++)
				exr -= s_ij * lmo_s[k][j] *
					fr_i->xr_fock_mat[fock_idx(i, k)];
			for (int l = 0; l < fr_j->n_lmo; l++)
				exr -= s_ij * lmo_s[i][l] *
					fr_j->xr_fock_mat[fock_idx(j, l)];
			exr += 2.0 * s_ij * t_ij;

			/* xr - third part */
			for (int jj = 0; jj < fr_j->n_atoms; jj++) {
				struct efp_atom *at = fr_j->atoms + jj;
				double r = vec_dist(fr_i->lmo_centroids + i,
							VEC(at->x));
				exr -= s_ij * s_ij * valence(at->znuc) / r;
			}
			for (int l = 0; l < fr_j->n_lmo; l++) {
				double r = vec_dist(fr_i->lmo_centroids + i,
						    fr_j->lmo_centroids + l);
				exr += 2.0 * s_ij * s_ij / r;
			}
			for (int ii = 0; ii < fr_i->n_atoms; ii++) {
				struct efp_atom *at = fr_i->atoms + ii;
				double r = vec_dist(fr_j->lmo_centroids + j,
							VEC(at->x));
				exr -= s_ij * s_ij * valence(at->znuc) / r;
			}
			for (int k = 0; k < fr_i->n_lmo; k++) {
				double r = vec_dist(fr_i->lmo_centroids + k,
						    fr_j->lmo_centroids + j);
				exr += 2.0 * s_ij * s_ij / r;
			}
			exr -= s_ij * s_ij / r_ij;
		}
	}

	*exr_out = 2.0 * exr;
	*ecp_out = ecp;

	if (efp->do_gradient)
		frag_frag_xr_grad(efp, frag_i, frag_j);
}

enum efp_result
efp_compute_xr(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_XR))
		return EFP_RESULT_SUCCESS;

	double exr = 0.0;
	double ecp = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:exr,ecp)
	for (int i = 0; i < efp->n_frag; i++) {
		for (int j = i + 1; j < efp->n_frag; j++) {
			double exr_out, ecp_out;

			frag_frag_xr(efp, i, j, &exr_out, &ecp_out);

			exr += exr_out;
			ecp += ecp_out;
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

	rotate_t2(rotmat, full_in, full_out);

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

	rotate_t3(rotmat, full_in, full_out);

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

void
efp_update_xr(struct frag *frag)
{
	const mat_t *rotmat = &frag->rotmat;

	/* update LMO centroids */
	for (int i = 0; i < frag->n_lmo; i++)
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			frag->lib->lmo_centroids + i, frag->lmo_centroids + i);

	/* update shells */
	for (int i = 0; i < frag->n_xr_shells; i++)
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			VEC(frag->lib->xr_shells[i].x),
			VEC(frag->xr_shells[i].x));

	/* rotate wavefunction */
	for (int lmo = 0; lmo < frag->n_lmo; lmo++) {
		const double *in = frag->lib->xr_wf + lmo * frag->xr_wf_size;
		double *out = frag->xr_wf + lmo * frag->xr_wf_size;

		for (int i = 0, func = 0; i < frag->n_xr_shells; i++) {
			switch (frag->xr_shells[i].type) {
			case 'S':
				func++;
				break;
			case 'L':
				func++;
				/* fall through */
			case 'P':
				mat_vec(rotmat, (const vec_t *)(in + func),
						(vec_t *)(out + func));
				func += 3;
				break;
			case 'D':
				rotate_func_d(rotmat, in + func, out + func);
				func += 6;
				break;
			case 'F':
				rotate_func_f(rotmat, in + func, out + func);
				func += 10;
				break;
			}
		}
	}
}
