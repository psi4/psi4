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

#include "balance.h"
#include "elec.h"
#include "private.h"

#define POL_SCF_TOL 1.0e-10
#define POL_SCF_MAX_ITER 80

double efp_get_pol_damp_tt(double, double, double);
enum efp_result efp_compute_id_direct(struct efp *);

struct id_work_data {
	double conv;
	vec_t *id_new;
	vec_t *id_conj_new;
};

double
efp_get_pol_damp_tt(double r, double pa, double pb)
{
	double ab = sqrt(pa * pb);
	double r2 = r * r;

	return (1.0 - exp(-ab * r2) * (1.0 + ab * r2));
}

static double
efp_get_pol_damp_tt_grad(double r, double pa, double pb)
{
	double ab = sqrt(pa * pb);
	double r2 = r * r;

	return (-2.0 * exp(-ab * r2) * (ab * ab * r2));
}

static vec_t
get_multipole_field(const vec_t *xyz, const struct multipole_pt *mult_pt,
			const struct swf *swf)
{
	vec_t field = vec_zero;

	vec_t dr = {
		xyz->x - mult_pt->x - swf->cell.x,
		xyz->y - mult_pt->y - swf->cell.y,
		xyz->z - mult_pt->z - swf->cell.z
	};

	double t1, t2;
	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;
	double r7 = r5 * r * r;

	/* charge */
	field.x += swf->swf * mult_pt->monopole * dr.x / r3;
	field.y += swf->swf * mult_pt->monopole * dr.y / r3;
	field.z += swf->swf * mult_pt->monopole * dr.z / r3;

	/* dipole */
	t1 = vec_dot(&mult_pt->dipole, &dr);

	field.x += swf->swf * (3.0 / r5 * t1 * dr.x - mult_pt->dipole.x / r3);
	field.y += swf->swf * (3.0 / r5 * t1 * dr.y - mult_pt->dipole.y / r3);
	field.z += swf->swf * (3.0 / r5 * t1 * dr.z - mult_pt->dipole.z / r3);

	/* quadrupole */
	t1 = quadrupole_sum(mult_pt->quadrupole, &dr);

	t2 = mult_pt->quadrupole[quad_idx(0, 0)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 0)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 0)] * dr.z;
	field.x += swf->swf * (-2.0 / r5 * t2 + 5.0 / r7 * t1 * dr.x);

	t2 = mult_pt->quadrupole[quad_idx(0, 1)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 1)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 1)] * dr.z;
	field.y += swf->swf * (-2.0 / r5 * t2 + 5.0 / r7 * t1 * dr.y);

	t2 = mult_pt->quadrupole[quad_idx(0, 2)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 2)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 2)] * dr.z;
	field.z += swf->swf * (-2.0 / r5 * t2 + 5.0 / r7 * t1 * dr.z);

	/* octupole-polarizability interactions are ignored */

	return (field);
}

static vec_t
get_elec_field(const struct efp *efp, size_t frag_idx, size_t pt_idx)
{
	const struct frag *fr_j = efp->frags + frag_idx;
	const struct polarizable_pt *pt = fr_j->polarizable_pts + pt_idx;
	vec_t elec_field = vec_zero;

	for (size_t i = 0; i < efp->n_frag; i++) {
		if (i == frag_idx || efp_skip_frag_pair(efp, i, frag_idx))
			continue;

		const struct frag *fr_i = efp->frags + i;
		struct swf swf = efp_make_swf(efp, fr_i, fr_j);

		/* field due to nuclei */
		for (size_t j = 0; j < fr_i->n_atoms; j++) {
			const struct efp_atom *at = fr_i->atoms + j;

			vec_t dr = {
				pt->x - at->x - swf.cell.x,
				pt->y - at->y - swf.cell.y,
				pt->z - at->z - swf.cell.z
			};

			double r = vec_len(&dr);
			double r3 = r * r * r;
			double p1 = 1.0;

			if (efp->opts.pol_damp == EFP_POL_DAMP_TT)
				p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp, fr_j->pol_damp);

			elec_field.x += swf.swf * at->znuc * dr.x / r3 * p1;
			elec_field.y += swf.swf * at->znuc * dr.y / r3 * p1;
			elec_field.z += swf.swf * at->znuc * dr.z / r3 * p1;
		}

		/* field due to multipoles */
		for (size_t j = 0; j < fr_i->n_multipole_pts; j++) {
			const struct multipole_pt *mult_pt = fr_i->multipole_pts + j;
			vec_t mult_field = get_multipole_field(CVEC(pt->x), mult_pt, &swf);

			vec_t dr = {
				pt->x - mult_pt->x - swf.cell.x,
				pt->y - mult_pt->y - swf.cell.y,
				pt->z - mult_pt->z - swf.cell.z
			};

			double r = vec_len(&dr);
			double p1 = 1.0;

			if (efp->opts.pol_damp == EFP_POL_DAMP_TT)
				p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp, fr_j->pol_damp);

			elec_field.x += mult_field.x * p1;
			elec_field.y += mult_field.y * p1;
			elec_field.z += mult_field.z * p1;
		}
	}

	if (efp->opts.terms & EFP_TERM_AI_POL) {
		/* field due to nuclei from ab initio subsystem */
		for (size_t i = 0; i < efp->n_ptc; i++) {
			vec_t dr = vec_sub(CVEC(pt->x), efp->ptc_xyz + i);

			double r = vec_len(&dr);
			double r3 = r * r * r;

			elec_field.x += efp->ptc[i] * dr.x / r3;
			elec_field.y += efp->ptc[i] * dr.y / r3;
			elec_field.z += efp->ptc[i] * dr.z / r3;
		}
	}

	return (elec_field);
}

static enum efp_result
add_electron_density_field(struct efp *efp)
{
	enum efp_result res;
	vec_t *xyz, *field;

	if (efp->get_electron_density_field == NULL)
		return (EFP_RESULT_SUCCESS);

	xyz = (vec_t *)malloc(efp->n_polarizable_pts * sizeof(vec_t));
	field = (vec_t *)malloc(efp->n_polarizable_pts * sizeof(vec_t));

	for (size_t i = 0, idx = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++, idx++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			xyz[idx].x = pt->x;
			xyz[idx].y = pt->y;
			xyz[idx].z = pt->z;
		}
	}

	if ((res = efp->get_electron_density_field(efp->n_polarizable_pts,
			(const double *)xyz, (double *)field,
			efp->get_electron_density_field_user_data)))
		goto error;

	for (size_t i = 0, idx = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++, idx++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			pt->elec_field_wf = field[idx];
		}
	}

error:
	free(xyz);
	free(field);
	return (res);
}

static void
compute_elec_field_range(struct efp *efp, size_t from, size_t to, void *data)
{
	vec_t *elec_field = (vec_t *)data;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (size_t i = from; i < to; i++) {
		const struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++)
			elec_field[frag->polarizable_offset + j] = get_elec_field(efp, i, j);
	}
}

static enum efp_result
compute_elec_field(struct efp *efp)
{
	vec_t *elec_field;
	enum efp_result res;

	elec_field = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));
	efp_balance_work(efp, compute_elec_field_range, elec_field);
	efp_allreduce((double *)elec_field, 3 * efp->n_polarizable_pts);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			pt->elec_field = elec_field[frag->polarizable_offset + j];
			pt->elec_field_wf = vec_zero;
		}
	}

	free(elec_field);

	if (efp->opts.terms & EFP_TERM_AI_POL)
		if ((res = add_electron_density_field(efp)))
			return (res);

	return (EFP_RESULT_SUCCESS);
}

static void
get_induced_dipole_field(struct efp *efp, size_t frag_idx,
			 struct polarizable_pt *pt,
			 vec_t *field, vec_t *field_conj)
{
	struct frag *fr_i = efp->frags + frag_idx;

	*field = vec_zero;
	*field_conj = vec_zero;

	for (size_t j = 0; j < efp->n_frag; j++) {
		if (j == frag_idx || efp_skip_frag_pair(efp, frag_idx, j))
			continue;

		struct frag *fr_j = efp->frags + j;
		struct swf swf = efp_make_swf(efp, fr_i, fr_j);

		for (size_t jj = 0; jj < fr_j->n_polarizable_pts; jj++) {
			struct polarizable_pt *pt_j = fr_j->polarizable_pts + jj;
			size_t idx = fr_j->polarizable_offset + jj;

			vec_t dr = {
				pt->x - pt_j->x + swf.cell.x,
				pt->y - pt_j->y + swf.cell.y,
				pt->z - pt_j->z + swf.cell.z
			};

			double r = vec_len(&dr);
			double r3 = r * r * r;
			double r5 = r3 * r * r;

			double t1 = vec_dot(&efp->indip[idx], &dr);
			double t2 = vec_dot(&efp->indipconj[idx], &dr);

			double p1 = 1.0;

			if (efp->opts.pol_damp == EFP_POL_DAMP_TT)
				p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp, fr_j->pol_damp);

			field->x -= swf.swf * p1 * (efp->indip[idx].x / r3 -
						3.0 * t1 * dr.x / r5);
			field->y -= swf.swf * p1 * (efp->indip[idx].y / r3 -
						3.0 * t1 * dr.y / r5);
			field->z -= swf.swf * p1 * (efp->indip[idx].z / r3 -
						3.0 * t1 * dr.z / r5);

			field_conj->x -= swf.swf * p1 * (efp->indipconj[idx].x / r3 -
						3.0 * t2 * dr.x / r5);
			field_conj->y -= swf.swf * p1 * (efp->indipconj[idx].y / r3 -
						3.0 * t2 * dr.y / r5);
			field_conj->z -= swf.swf * p1 * (efp->indipconj[idx].z / r3 -
						3.0 * t2 * dr.z / r5);
		}
	}
}

static void
compute_id_range(struct efp *efp, size_t from, size_t to, void *data)
{
	double conv = 0.0;
	vec_t *id_new, *id_conj_new;

	id_new = ((struct id_work_data *)data)->id_new;
	id_conj_new = ((struct id_work_data *)data)->id_conj_new;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:conv)
#endif
	for (size_t i = from; i < to; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			size_t idx = frag->polarizable_offset + j;
			vec_t field, field_conj;

			/* electric field from other induced dipoles */
			get_induced_dipole_field(efp, i, pt, &field, &field_conj);

			/* add field that doesn't change during scf */
			field.x += pt->elec_field.x + pt->elec_field_wf.x;
			field.y += pt->elec_field.y + pt->elec_field_wf.y;
			field.z += pt->elec_field.z + pt->elec_field_wf.z;

			field_conj.x += pt->elec_field.x + pt->elec_field_wf.x;
			field_conj.y += pt->elec_field.y + pt->elec_field_wf.y;
			field_conj.z += pt->elec_field.z + pt->elec_field_wf.z;

			id_new[idx] = mat_vec(&pt->tensor, &field);
			id_conj_new[idx] = mat_trans_vec(&pt->tensor, &field_conj);

			conv += vec_dist(&id_new[idx], &efp->indip[idx]);
			conv += vec_dist(&id_conj_new[idx], &efp->indipconj[idx]);
		}
	}

	((struct id_work_data *)data)->conv += conv;
}

static double
pol_scf_iter(struct efp *efp)
{
	struct id_work_data data;

	data.conv = 0.0;
	data.id_new = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));
	data.id_conj_new = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));

	efp_balance_work(efp, compute_id_range, &data);

	efp_allreduce((double *)data.id_new, 3 * efp->n_polarizable_pts);
	efp_allreduce((double *)data.id_conj_new, 3 * efp->n_polarizable_pts);
	efp_allreduce(&data.conv, 1);

	memcpy(efp->indip, data.id_new, efp->n_polarizable_pts * sizeof(vec_t));
	memcpy(efp->indipconj, data.id_conj_new, efp->n_polarizable_pts * sizeof(vec_t));

	free(data.id_new);
	free(data.id_conj_new);

	return (data.conv / efp->n_polarizable_pts / 2);
}

static void
compute_energy_range(struct efp *efp, size_t from, size_t to, void *data)
{
	double energy = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:energy)
#endif
	for (size_t i = from; i < to; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			size_t idx = frag->polarizable_offset + j;

			energy += 0.5 * vec_dot(&efp->indipconj[idx],
						&pt->elec_field_wf) -
				  0.5 * vec_dot(&efp->indip[idx],
						&pt->elec_field);
		}
	}

	*(double *)data += energy;
}

static enum efp_result
efp_compute_id_iterative(struct efp *efp)
{
	memset(efp->indip, 0, efp->n_polarizable_pts * sizeof(vec_t));
	memset(efp->indipconj, 0, efp->n_polarizable_pts * sizeof(vec_t));

	for (size_t iter = 1; iter <= POL_SCF_MAX_ITER; iter++) {
		if (pol_scf_iter(efp) < POL_SCF_TOL)
			break;

		if (iter == POL_SCF_MAX_ITER)
			return (EFP_RESULT_POL_NOT_CONVERGED);
	}

	return (EFP_RESULT_SUCCESS);
}

enum efp_result
efp_compute_pol_energy(struct efp *efp, double *energy)
{
	enum efp_result res;

	assert(energy);

	if ((res = compute_elec_field(efp)))
		return (res);

	switch (efp->opts.pol_driver) {
		case EFP_POL_DRIVER_ITERATIVE:
			res = efp_compute_id_iterative(efp);
			break;
		case EFP_POL_DRIVER_DIRECT:
			res = efp_compute_id_direct(efp);
			break;
	}

	if (res)
		return (res);

	*energy = 0.0;
	efp_balance_work(efp, compute_energy_range, energy);
	efp_allreduce(energy, 1);

	return (EFP_RESULT_SUCCESS);
}

static void
compute_grad_point(struct efp *efp, size_t frag_idx, size_t pt_idx)
{
	const struct frag *fr_i = efp->frags + frag_idx;
	const struct polarizable_pt *pt_i = fr_i->polarizable_pts + pt_idx;
	size_t idx_i = fr_i->polarizable_offset + pt_idx;

	vec_t dipole_i = {
		0.5 * (efp->indip[idx_i].x + efp->indipconj[idx_i].x),
		0.5 * (efp->indip[idx_i].y + efp->indipconj[idx_i].y),
		0.5 * (efp->indip[idx_i].z + efp->indipconj[idx_i].z)
	};

	for (size_t j = 0; j < efp->n_frag; j++) {
		if (j == frag_idx || efp_skip_frag_pair(efp, frag_idx, j))
			continue;

		struct frag *fr_j = efp->frags + j;
		struct swf swf = efp_make_swf(efp, fr_i, fr_j);

		/* energy without switching applied */
		double energy = 0.0;

		/* induced dipole - nuclei */
		for (size_t k = 0; k < fr_j->n_atoms; k++) {
			struct efp_atom *at_j = fr_j->atoms + k;

			vec_t dr = {
				at_j->x - pt_i->x - swf.cell.x,
				at_j->y - pt_i->y - swf.cell.y,
				at_j->z - pt_i->z - swf.cell.z
			};

			double p1 = 1.0, p2 = 0.0;

			if (efp->opts.pol_damp == EFP_POL_DAMP_TT) {
				double r = vec_len(&dr);

				p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp,
						fr_j->pol_damp);
				p2 = efp_get_pol_damp_tt_grad(r, fr_i->pol_damp,
						fr_j->pol_damp);
			}

			vec_t force, add_i, add_j;

			double e = -efp_charge_dipole_energy(at_j->znuc, &dipole_i, &dr);

			efp_charge_dipole_grad(at_j->znuc, &dipole_i, &dr,
					       &force, &add_j, &add_i);
			vec_negate(&force);

			vec_scale(&force, p1);
			vec_scale(&add_i, p1);
			vec_scale(&add_j, p1);

			force.x += p2 * e * dr.x;
			force.y += p2 * e * dr.y;
			force.z += p2 * e * dr.z;

			vec_scale(&force, swf.swf);
			vec_scale(&add_i, swf.swf);
			vec_scale(&add_j, swf.swf);

			efp_add_force(efp->grad + frag_idx, CVEC(fr_i->x),
					CVEC(pt_i->x), &force, &add_i);
			efp_sub_force(efp->grad + j, CVEC(fr_j->x),
					CVEC(at_j->x), &force, &add_j);
			efp_add_stress(&swf.dr, &force, &efp->stress);

			energy += p1 * e;
		}

		/* induced dipole - multipoles */
		for (size_t k = 0; k < fr_j->n_multipole_pts; k++) {
			struct multipole_pt *pt_j = fr_j->multipole_pts + k;

			vec_t dr = {
				pt_j->x - pt_i->x - swf.cell.x,
				pt_j->y - pt_i->y - swf.cell.y,
				pt_j->z - pt_i->z - swf.cell.z
			};

			double p1 = 1.0, p2 = 0.0;

			if (efp->opts.pol_damp == EFP_POL_DAMP_TT) {
				double r = vec_len(&dr);

				p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp,
						fr_j->pol_damp);
				p2 = efp_get_pol_damp_tt_grad(r, fr_i->pol_damp,
						fr_j->pol_damp);
			}

			double e = 0.0;

			vec_t force_, add_i_, add_j_;
			vec_t force = vec_zero, add_i = vec_zero, add_j = vec_zero;

			/* induced dipole - charge */
			e -= efp_charge_dipole_energy(pt_j->monopole, &dipole_i, &dr);

			efp_charge_dipole_grad(pt_j->monopole, &dipole_i, &dr,
					       &force_, &add_j_, &add_i_);
			vec_negate(&force_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			/* induced dipole - dipole */
			e += efp_dipole_dipole_energy(&dipole_i, &pt_j->dipole, &dr);

			efp_dipole_dipole_grad(&dipole_i, &pt_j->dipole, &dr,
					       &force_, &add_i_, &add_j_);
			vec_negate(&add_j_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			/* induced dipole - quadrupole */
			e += efp_dipole_quadrupole_energy(&dipole_i, pt_j->quadrupole, &dr);

			efp_dipole_quadrupole_grad(&dipole_i, pt_j->quadrupole, &dr,
						   &force_, &add_i_, &add_j_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			/* induced dipole - octupole interactions are ignored */

			vec_scale(&force, p1);
			vec_scale(&add_i, p1);
			vec_scale(&add_j, p1);

			force.x += p2 * e * dr.x;
			force.y += p2 * e * dr.y;
			force.z += p2 * e * dr.z;

			vec_scale(&force, swf.swf);
			vec_scale(&add_i, swf.swf);
			vec_scale(&add_j, swf.swf);

			efp_add_force(efp->grad + frag_idx, CVEC(fr_i->x),
					CVEC(pt_i->x), &force, &add_i);
			efp_sub_force(efp->grad + j, CVEC(fr_j->x),
					CVEC(pt_j->x), &force, &add_j);
			efp_add_stress(&swf.dr, &force, &efp->stress);

			energy += p1 * e;
		}

		/* induced dipole - induced dipoles */
		for (size_t jj = 0; jj < fr_j->n_polarizable_pts; jj++) {
			struct polarizable_pt *pt_j = fr_j->polarizable_pts + jj;
			size_t idx_j = fr_j->polarizable_offset + jj;

			vec_t dr = {
				pt_j->x - pt_i->x - swf.cell.x,
				pt_j->y - pt_i->y - swf.cell.y,
				pt_j->z - pt_i->z - swf.cell.z
			};

			vec_t half_dipole_i = {
				0.5 * efp->indip[idx_i].x,
				0.5 * efp->indip[idx_i].y,
				0.5 * efp->indip[idx_i].z
			};

			double p1 = 1.0, p2 = 0.0;

			if (efp->opts.pol_damp == EFP_POL_DAMP_TT) {
				double r = vec_len(&dr);

				p1 = efp_get_pol_damp_tt(r, fr_i->pol_damp,
						fr_j->pol_damp);
				p2 = efp_get_pol_damp_tt_grad(r, fr_i->pol_damp,
						fr_j->pol_damp);
			}

			vec_t force, add_i, add_j;

			double e = efp_dipole_dipole_energy(&half_dipole_i,
						&efp->indipconj[idx_j], &dr);

			efp_dipole_dipole_grad(&half_dipole_i, &efp->indipconj[idx_j],
						&dr, &force, &add_i, &add_j);
			vec_negate(&add_j);

			vec_scale(&force, p1);
			vec_scale(&add_i, p1);
			vec_scale(&add_j, p1);

			force.x += p2 * e * dr.x;
			force.y += p2 * e * dr.y;
			force.z += p2 * e * dr.z;

			vec_scale(&force, swf.swf);
			vec_scale(&add_i, swf.swf);
			vec_scale(&add_j, swf.swf);

			efp_add_force(efp->grad + frag_idx, CVEC(fr_i->x),
					CVEC(pt_i->x), &force, &add_i);
			efp_sub_force(efp->grad + j, CVEC(fr_j->x),
					CVEC(pt_j->x), &force, &add_j);
			efp_add_stress(&swf.dr, &force, &efp->stress);

			energy += p1 * e;
		}

		vec_t force = {
			swf.dswf.x * energy,
			swf.dswf.y * energy,
			swf.dswf.z * energy
		};

		six_atomic_add_xyz(efp->grad + frag_idx, &force);
		six_atomic_sub_xyz(efp->grad + j, &force);
		efp_add_stress(&swf.dr, &force, &efp->stress);
	}

	/* induced dipole - ab initio nuclei */
	if (efp->opts.terms & EFP_TERM_AI_POL) {
		for (size_t j = 0; j < efp->n_ptc; j++) {
			vec_t dr = vec_sub(efp->ptc_xyz + j, CVEC(pt_i->x));
			vec_t force, add_i, add_j;

			efp_charge_dipole_grad(efp->ptc[j], &dipole_i, &dr,
					       &force, &add_j, &add_i);
			vec_negate(&add_i);

			vec_atomic_add(efp->ptc_grad + j, &force);
			efp_sub_force(efp->grad + frag_idx, CVEC(fr_i->x),
					CVEC(pt_i->x), &force, &add_i);
		}
	}
}

static void
compute_grad_range(struct efp *efp, size_t from, size_t to, void *data)
{
	(void)data;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (size_t i = from; i < to; i++)
		for (size_t j = 0; j < efp->frags[i].n_polarizable_pts; j++)
			compute_grad_point(efp, i, j);
}

enum efp_result
efp_compute_pol(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_POL) && !(efp->opts.terms & EFP_TERM_AI_POL))
		return EFP_RESULT_SUCCESS;

	enum efp_result res;

	if ((res = efp_compute_pol_energy(efp, &efp->energy.polarization)))
		return res;

	if (efp->do_gradient)
		efp_balance_work(efp, compute_grad_range, NULL);

	return EFP_RESULT_SUCCESS;
}

void
efp_update_pol(struct frag *frag)
{
	for (size_t i = 0; i < frag->n_polarizable_pts; i++) {
		efp_move_pt(CVEC(frag->x), &frag->rotmat,
			CVEC(frag->lib->polarizable_pts[i].x),
			VEC(frag->polarizable_pts[i].x));

		const mat_t *in = &frag->lib->polarizable_pts[i].tensor;
		mat_t *out = &frag->polarizable_pts[i].tensor;

		efp_rotate_t2(&frag->rotmat, (const double *)in, (double *)out);
	}
}

EFP_EXPORT enum efp_result
efp_get_electric_field(struct efp *efp, size_t frag_idx, const double *xyz, double *field)
{
	assert(efp);
	assert(frag_idx < efp->n_frag);
	assert(xyz);
	assert(field);

	const struct frag *frag = efp->frags + frag_idx;
	vec_t elec_field = vec_zero;

	for (size_t i = 0; i < efp->n_frag; i++) {
		if (i == frag_idx || efp_skip_frag_pair(efp, i, frag_idx))
			continue;

		const struct frag *fr_i = efp->frags + i;
		struct swf swf = efp_make_swf(efp, fr_i, frag);

		/* field due to nuclei */
		for (size_t j = 0; j < fr_i->n_atoms; j++) {
			const struct efp_atom *at = fr_i->atoms + j;

			vec_t dr = {
				xyz[0] - at->x - swf.cell.x,
				xyz[1] - at->y - swf.cell.y,
				xyz[2] - at->z - swf.cell.z
			};

			double r = vec_len(&dr);
			double r3 = r * r * r;

			elec_field.x += swf.swf * at->znuc * dr.x / r3;
			elec_field.y += swf.swf * at->znuc * dr.y / r3;
			elec_field.z += swf.swf * at->znuc * dr.z / r3;
		}

		/* field due to multipoles */
		for (size_t j = 0; j < fr_i->n_multipole_pts; j++) {
			const struct multipole_pt *mpt = fr_i->multipole_pts + j;
			vec_t mult_field = get_multipole_field((const vec_t *)xyz, mpt, &swf);

			elec_field.x += mult_field.x;
			elec_field.y += mult_field.y;
			elec_field.z += mult_field.z;
		}

		/* field due to induced dipoles */
		for (size_t j = 0; j < fr_i->n_polarizable_pts; j++) {
			struct polarizable_pt *pt_i = fr_i->polarizable_pts + j;
			size_t idx = fr_i->polarizable_offset + j;

			vec_t dr = {
				xyz[0] - pt_i->x - swf.cell.x,
				xyz[1] - pt_i->y - swf.cell.y,
				xyz[2] - pt_i->z - swf.cell.z
			};

			double r = vec_len(&dr);
			double r3 = r * r * r;
			double r5 = r3 * r * r;
			double t1 = vec_dot(&efp->indip[idx], &dr);

			elec_field.x -= swf.swf * (efp->indip[idx].x / r3 -
						3.0 * t1 * dr.x / r5);
			elec_field.y -= swf.swf * (efp->indip[idx].y / r3 -
						3.0 * t1 * dr.y / r5);
			elec_field.z -= swf.swf * (efp->indip[idx].z / r3 -
						3.0 * t1 * dr.z / r5);
		}
	}

	if (efp->opts.terms & EFP_TERM_AI_POL) {
		/* field due to nuclei from ab initio subsystem */
		for (size_t i = 0; i < efp->n_ptc; i++) {
			vec_t dr = vec_sub((const vec_t *)xyz, efp->ptc_xyz + i);

			double r = vec_len(&dr);
			double r3 = r * r * r;

			elec_field.x += efp->ptc[i] * dr.x / r3;
			elec_field.y += efp->ptc[i] * dr.y / r3;
			elec_field.z += efp->ptc[i] * dr.z / r3;
		}
	}

	*((vec_t *)field) = elec_field;
	return (EFP_RESULT_SUCCESS);
}
