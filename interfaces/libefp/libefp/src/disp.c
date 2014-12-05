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

#include "private.h"

static const double weights[] = {
	0.72086099022968040154E-02, 0.17697067815034886394E-01,
	0.30660908596251749739E-01, 0.48381293256249884995E-01,
	0.74878830420650517080E-01, 0.11806515901361630228E+00,
	0.19535413832209084204E+00, 0.35055692324483221824E+00,
	0.71577113554429568336E+00, 0.18140975997632396972E+01,
	0.69792344511487082324E+01, 0.83248093882965845391E+02
};

static double
get_damp_tt(double r)
{
	static const double a = 1.5; /* Tang-Toennies damping parameter */

	double ra = r * a;
	double ra2 = ra * ra;
	double ra3 = ra2 * ra;
	double ra4 = ra3 * ra;
	double ra5 = ra4 * ra;
	double ra6 = ra5 * ra;

	return 1.0 - exp(-ra) * (1.0 + ra + ra2 / 2.0 + ra3 / 6.0 +
			ra4 / 24.0 + ra5 / 120.0 + ra6 / 720.0);
}

static double
get_damp_tt_grad(double r)
{
	static const double a = 1.5; /* Tang-Toennies damping parameter */

	double ra = r * a;
	double ra2 = ra * ra;
	double ra6 = ra2 * ra2 * ra2;

	return a * exp(-ra) * ra6 / 720.0;
}

static double
disp_tt(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
	size_t pt_i_idx, size_t pt_j_idx, double sum, const struct swf *swf)
{
	const struct frag *fr_i = efp->frags + fr_i_idx;
	const struct frag *fr_j = efp->frags + fr_j_idx;

	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - pt_i->x - swf->cell.x,
		pt_j->y - pt_i->y - swf->cell.y,
		pt_j->z - pt_i->z - swf->cell.z
	};

	double r = vec_len(&dr);
	double r2 = r * r;
	double r6 = r2 * r2 * r2;

	double damp = get_damp_tt(r);
	double energy = -4.0 / 3.0 * sum * damp / r6;

	if (efp->do_gradient) {
		double gdamp = get_damp_tt_grad(r);
		double g = 4.0 / 3.0 * sum * (gdamp / r - 6.0 * damp / r2) / r6;

		vec_t force = {
			g * dr.x * swf->swf,
			g * dr.y * swf->swf,
			g * dr.z * swf->swf
		};

		efp_add_force(efp->grad + fr_i_idx, CVEC(fr_i->x),
				CVEC(pt_i->x), &force, NULL);
		efp_sub_force(efp->grad + fr_j_idx, CVEC(fr_j->x),
				CVEC(pt_j->x), &force, NULL);
		efp_add_stress(&swf->dr, &force, &efp->stress);
	}

	return energy;
}

static double
disp_overlap(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
	     size_t pt_i_idx, size_t pt_j_idx, double s_ij, six_t ds_ij,
	     double sum, const struct swf *swf)
{
	const struct frag *fr_i = efp->frags + fr_i_idx;
	const struct frag *fr_j = efp->frags + fr_j_idx;

	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - pt_i->x - swf->cell.x,
		pt_j->y - pt_i->y - swf->cell.y,
		pt_j->z - pt_i->z - swf->cell.z
	};

	double r = vec_len(&dr);
	double r2 = r * r;
	double r6 = r2 * r2 * r2;

	double ln_s = 0.0;
	double damp = 1.0;

	if (fabs(s_ij) > 1.0e-5) {
		ln_s = log(fabs(s_ij));
		damp = 1.0 - s_ij * s_ij * (1.0 - 2.0 * ln_s + 2.0 * ln_s * ln_s);
	}

	double energy = -4.0 / 3.0 * sum * damp / r6;

	if (efp->do_gradient) {
		vec_t force, torque_i, torque_j;

		double t1 = -8.0 * sum / r6 / r2 * damp;
		double t2 = -16.0 / 3.0 * sum / r6 * ln_s * ln_s * s_ij;

		vec_t dr_i = vec_sub(CVEC(pt_i->x), VEC(fr_i->x));
		vec_t dr_j = vec_sub(CVEC(pt_j->x), VEC(fr_j->x));

		force.x = (t1 * dr.x - t2 * ds_ij.x) * swf->swf;
		force.y = (t1 * dr.y - t2 * ds_ij.y) * swf->swf;
		force.z = (t1 * dr.z - t2 * ds_ij.z) * swf->swf;

		torque_i.x = swf->swf * (t1 * (dr.z * dr_i.y - dr.y * dr_i.z) +
						t2 * ds_ij.a);
		torque_i.y = swf->swf * (t1 * (dr.x * dr_i.z - dr.z * dr_i.x) +
						t2 * ds_ij.b);
		torque_i.z = swf->swf * (t1 * (dr.y * dr_i.x - dr.x * dr_i.y) +
						t2 * ds_ij.c);

		torque_j.x = swf->swf * (t1 * (dr.z * dr_j.y - dr.y * dr_j.z) +
			     t2 * (ds_ij.z * swf->dr.y - ds_ij.y * swf->dr.z) +
			     t2 * ds_ij.a);
		torque_j.y = swf->swf * (t1 * (dr.x * dr_j.z - dr.z * dr_j.x) +
			     t2 * (ds_ij.x * swf->dr.z - ds_ij.z * swf->dr.x) +
			     t2 * ds_ij.b);
		torque_j.z = swf->swf * (t1 * (dr.y * dr_j.x - dr.x * dr_j.y) +
			     t2 * (ds_ij.y * swf->dr.x - ds_ij.x * swf->dr.y) +
			     t2 * ds_ij.c);

		six_atomic_add_xyz(efp->grad + fr_i_idx, &force);
		six_atomic_add_abc(efp->grad + fr_i_idx, &torque_i);

		six_atomic_sub_xyz(efp->grad + fr_j_idx, &force);
		six_atomic_sub_abc(efp->grad + fr_j_idx, &torque_j);

		efp_add_stress(&swf->dr, &force, &efp->stress);
	}

	return energy;
}

static double
disp_off(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
	 size_t pt_i_idx, size_t pt_j_idx, double sum, const struct swf *swf)
{
	const struct frag *fr_i = efp->frags + fr_i_idx;
	const struct frag *fr_j = efp->frags + fr_j_idx;

	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - pt_i->x - swf->cell.x,
		pt_j->y - pt_i->y - swf->cell.y,
		pt_j->z - pt_i->z - swf->cell.z
	};

	double r = vec_len(&dr);
	double r2 = r * r;
	double r6 = r2 * r2 * r2;

	double energy = -4.0 / 3.0 * sum / r6;

	if (efp->do_gradient) {
		double r8 = r6 * r2;
		double g = -8.0 * sum / r8;

		vec_t force = {
			g * dr.x * swf->swf,
			g * dr.y * swf->swf,
			g * dr.z * swf->swf
		};

		efp_add_force(efp->grad + fr_i_idx, CVEC(fr_i->x),
				CVEC(pt_i->x), &force, NULL);
		efp_sub_force(efp->grad + fr_j_idx, CVEC(fr_j->x),
				CVEC(pt_j->x), &force, NULL);
		efp_add_stress(&swf->dr, &force, &efp->stress);
	}

	return energy;
}

static double
point_point_disp(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
		 size_t pt_i_idx, size_t pt_j_idx, double s, six_t ds,
		 const struct swf *swf)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	double sum = 0.0;

	for (size_t k = 0; k < ARRAY_SIZE(weights); k++) {
		double tr_i = (pt_i->tensor[k].xx + pt_i->tensor[k].yy + pt_i->tensor[k].zz) / 3;
		double tr_j = (pt_j->tensor[k].xx + pt_j->tensor[k].yy + pt_j->tensor[k].zz) / 3;

		sum += weights[k] * tr_i * tr_j;
	}

	switch (efp->opts.disp_damp) {
		case EFP_DISP_DAMP_TT:
			return disp_tt(efp, fr_i_idx, fr_j_idx,
				pt_i_idx, pt_j_idx, sum, swf);
		case EFP_DISP_DAMP_OVERLAP:
			return disp_overlap(efp, fr_i_idx, fr_j_idx,
				pt_i_idx, pt_j_idx, s, ds, sum, swf);
		case EFP_DISP_DAMP_OFF:
			return disp_off(efp, fr_i_idx, fr_j_idx,
				pt_i_idx, pt_j_idx, sum, swf);
	}
	assert(0);
}

/*
 * Reference:
 *
 * Ivana Adamovic, Mark Gordon
 *
 * Dynamic polarizability, dispersion coefficient C6 and dispersion energy in
 * the effective fragment potential method
 *
 * Mol. Phys. 103, 379 (2005)
 */
double
efp_frag_frag_disp(struct efp *efp, size_t frag_i, size_t frag_j,
			const double *s, const six_t *ds)
{
	double energy = 0.0;

	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	size_t n_disp_i = fr_i->n_dynamic_polarizable_pts;
	size_t n_disp_j = fr_j->n_dynamic_polarizable_pts;

	struct swf swf = efp_make_swf(efp, fr_i, fr_j);

	for (size_t ii = 0, idx = 0; ii < n_disp_i; ii++)
		for (size_t jj = 0; jj < n_disp_j; jj++, idx++)
			energy += point_point_disp(efp, frag_i, frag_j, ii, jj,
						   s[idx], ds[idx], &swf);

	vec_t force = {
		swf.dswf.x * energy,
		swf.dswf.y * energy,
		swf.dswf.z * energy
	};

	six_atomic_add_xyz(efp->grad + frag_i, &force);
	six_atomic_sub_xyz(efp->grad + frag_j, &force);
	efp_add_stress(&swf.dr, &force, &efp->stress);

	return energy * swf.swf;
}

void
efp_update_disp(struct frag *frag)
{
	for (size_t i = 0; i < frag->n_dynamic_polarizable_pts; i++) {
		const struct dynamic_polarizable_pt *pt_in =
					frag->lib->dynamic_polarizable_pts + i;
		struct dynamic_polarizable_pt *pt_out =
					frag->dynamic_polarizable_pts + i;

		efp_move_pt(CVEC(frag->x), &frag->rotmat, CVEC(pt_in->x), VEC(pt_out->x));

		for (size_t j = 0; j < 12; j++) {
			const mat_t *in = pt_in->tensor + j;
			mat_t *out = pt_out->tensor + j;

			efp_rotate_t2(&frag->rotmat, (const double *)in, (double *)out);
		}
	}
}
