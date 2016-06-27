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

#include "balance.h"
#include "elec.h"
#include "private.h"

static double
get_screen_damping(double r_ij, double pi, double pj)
{
	if (pj == HUGE_VAL) {   /* j is nucleus */
		return 1.0 - exp(-pi * r_ij);
	}
	else if (fabs(pi - pj) < 1.0e-5) {
		return 1.0 - (1.0 + 0.5 * pi * r_ij) * exp(-pi * r_ij);
	}
	else {
		return 1.0 - exp(-pi * r_ij) * pj * pj / (pj * pj - pi * pi) -
			     exp(-pj * r_ij) * pi * pi / (pi * pi - pj * pj);
	}
}

static double
get_screen_damping_grad(double r_ij, double pi, double pj)
{
	if (pj == HUGE_VAL) {   /* j is nucleus */
		return 1.0 - exp(-r_ij * pi) * (1.0 + pi * r_ij);
	}
	else if (fabs(pi - pj) < 1.0e-5) {
		return 1.0 - exp(-r_ij * pi) * (1.0 + pi * r_ij +
						0.5 * pi * pi * r_ij * r_ij);
	}
	else {
		return 1.0 - exp(-r_ij * pi) * (1.0 + pi * r_ij) *
						pj * pj / (pj * pj - pi * pi) -
			     exp(-r_ij * pj) * (1.0 + pj * r_ij) *
						pi * pi / (pi * pi - pj * pj);
	}
}

static double
atom_mult_energy(struct efp *efp, struct frag *fr_i, struct frag *fr_j,
		 size_t atom_i_idx, size_t pt_j_idx, const struct swf *swf)
{
	struct efp_atom *at_i = fr_i->atoms + atom_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - at_i->x - swf->cell.x,
		pt_j->y - at_i->y - swf->cell.y,
		pt_j->z - at_i->z - swf->cell.z
	};

	double energy = 0.0, ccdamp = 1.0;

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double r = vec_len(&dr);
		double sp = fr_j->screen_params[pt_j_idx];
		ccdamp = get_screen_damping(r, sp, HUGE_VAL);
	}

	/* charge - monopole */
	energy += ccdamp * efp_charge_charge_energy(at_i->znuc, pt_j->monopole, &dr);

	/* charge - dipole */
	energy += efp_charge_dipole_energy(at_i->znuc, &pt_j->dipole, &dr);

	/* charge - quadrupole */
	energy += efp_charge_quadrupole_energy(at_i->znuc, pt_j->quadrupole, &dr);

	/* charge - octupole */
	energy += efp_charge_octupole_energy(at_i->znuc, pt_j->octupole, &dr);

	return energy;
}

static void
atom_mult_grad(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
	       size_t atom_i_idx, size_t pt_j_idx, const struct swf *swf)
{
	const struct frag *fr_i = efp->frags + fr_i_idx;
	const struct frag *fr_j = efp->frags + fr_j_idx;
	const struct efp_atom *at_i = fr_i->atoms + atom_i_idx;
	const struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - at_i->x - swf->cell.x,
		pt_j->y - at_i->y - swf->cell.y,
		pt_j->z - at_i->z - swf->cell.z
	};

	vec_t force_, torque_i_, torque_j_;
	vec_t force = vec_zero, torque_i = vec_zero, torque_j = vec_zero;

	/* charge - charge */
	efp_charge_charge_grad(at_i->znuc, pt_j->monopole, &dr,
			       &force_, &torque_i_, &torque_j_);

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double r = vec_len(&dr);
		double sp = fr_j->screen_params[pt_j_idx];
		double gdamp = get_screen_damping_grad(r, sp, HUGE_VAL);

		force_.x *= gdamp;
		force_.y *= gdamp;
		force_.z *= gdamp;
	}

	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* charge - dipole */
	efp_charge_dipole_grad(at_i->znuc, &pt_j->dipole, &dr,
			       &force_, &torque_i_, &torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* charge - quadrupole */
	efp_charge_quadrupole_grad(at_i->znuc, pt_j->quadrupole, &dr,
				   &force_, &torque_i_, &torque_j_);
	vec_negate(&torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* charge - octupole */
	efp_charge_octupole_grad(at_i->znuc, pt_j->octupole, &dr,
				 &force_, &torque_i_, &torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	vec_scale(&force, swf->swf);
	vec_scale(&torque_i, swf->swf);
	vec_scale(&torque_j, swf->swf);

	efp_add_force(efp->grad + fr_i_idx, CVEC(fr_i->x), CVEC(at_i->x), &force, &torque_i);
	efp_sub_force(efp->grad + fr_j_idx, CVEC(fr_j->x), CVEC(pt_j->x), &force, &torque_j);
	efp_add_stress(&swf->dr, &force, &efp->stress);
}

static double
mult_mult_energy(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
		 size_t pt_i_idx, size_t pt_j_idx, const struct swf *swf)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - pt_i->x - swf->cell.x,
		pt_j->y - pt_i->y - swf->cell.y,
		pt_j->z - pt_i->z - swf->cell.z
	};

	double energy = 0.0, ccdamp = 1.0;

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double r = vec_len(&dr);
		double screen_i = fr_i->screen_params[pt_i_idx];
		double screen_j = fr_j->screen_params[pt_j_idx];
		ccdamp = get_screen_damping(r, screen_i, screen_j);
	}

	/* monopole - monopole */
	energy += ccdamp * efp_charge_charge_energy(pt_i->monopole, pt_j->monopole, &dr);

	/* monopole - dipole */
	energy += efp_charge_dipole_energy(pt_i->monopole, &pt_j->dipole, &dr);

	/* dipole - monopole */
	energy -= efp_charge_dipole_energy(pt_j->monopole, &pt_i->dipole, &dr);

	/* monopole - quadrupole */
	energy += efp_charge_quadrupole_energy(pt_i->monopole, pt_j->quadrupole, &dr);

	/* quadrupole - monopole */
	energy += efp_charge_quadrupole_energy(pt_j->monopole, pt_i->quadrupole, &dr);

	/* monopole - octupole */
	energy += efp_charge_octupole_energy(pt_i->monopole, pt_j->octupole, &dr);

	/* octupole - monopole */
	energy -= efp_charge_octupole_energy(pt_j->monopole, pt_i->octupole, &dr);

	/* dipole - dipole */
	energy += efp_dipole_dipole_energy(&pt_i->dipole, &pt_j->dipole, &dr);

	/* dipole - quadrupole */
	energy += efp_dipole_quadrupole_energy(&pt_i->dipole, pt_j->quadrupole, &dr);

	/* quadrupole - dipole */
	energy -= efp_dipole_quadrupole_energy(&pt_j->dipole, pt_i->quadrupole, &dr);

	/* quadrupole - quadrupole */
	energy += efp_quadrupole_quadrupole_energy(pt_i->quadrupole, pt_j->quadrupole, &dr);

	return energy;
}

static void
mult_mult_grad(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx,
	       size_t pt_i_idx, size_t pt_j_idx, const struct swf *swf)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = {
		pt_j->x - pt_i->x - swf->cell.x,
		pt_j->y - pt_i->y - swf->cell.y,
		pt_j->z - pt_i->z - swf->cell.z
	};

	vec_t force_, torque_i_, torque_j_;
	vec_t force = vec_zero, torque_i = vec_zero, torque_j = vec_zero;

	/* monopole - monopole */
	efp_charge_charge_grad(pt_i->monopole, pt_j->monopole, &dr,
			       &force_, &torque_i_, &torque_j_);

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double r = vec_len(&dr);

		double screen_i = fr_i->screen_params[pt_i_idx];
		double screen_j = fr_j->screen_params[pt_j_idx];

		double gdamp = get_screen_damping_grad(r, screen_i, screen_j);

		force_.x *= gdamp;
		force_.y *= gdamp;
		force_.z *= gdamp;
	}

	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* monopole - dipole */
	efp_charge_dipole_grad(pt_i->monopole, &pt_j->dipole, &dr,
			       &force_, &torque_i_, &torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* dipole - monopole */
	efp_charge_dipole_grad(pt_j->monopole, &pt_i->dipole, &dr,
			       &force_, &torque_j_, &torque_i_);
	vec_negate(&force_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* monopole - quadrupole */
	efp_charge_quadrupole_grad(pt_i->monopole, pt_j->quadrupole, &dr,
				   &force_, &torque_i_, &torque_j_);
	vec_negate(&torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* quadrupole - monopole */
	efp_charge_quadrupole_grad(pt_j->monopole, pt_i->quadrupole, &dr,
				   &force_, &torque_j_, &torque_i_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* monopole - octupole */
	efp_charge_octupole_grad(pt_i->monopole, pt_j->octupole, &dr,
				 &force_, &torque_i_, &torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* octupole - monopole */
	efp_charge_octupole_grad(pt_j->monopole, pt_i->octupole, &dr,
				 &force_, &torque_j_, &torque_i_);
	vec_negate(&force_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* dipole - dipole */
	efp_dipole_dipole_grad(&pt_i->dipole, &pt_j->dipole, &dr,
			       &force_, &torque_i_, &torque_j_);
	vec_negate(&torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* dipole - quadrupole */
	efp_dipole_quadrupole_grad(&pt_i->dipole, pt_j->quadrupole, &dr,
				   &force_, &torque_i_, &torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* quadrupole - dipole */
	efp_dipole_quadrupole_grad(&pt_j->dipole, pt_i->quadrupole, &dr,
				   &force_, &torque_j_, &torque_i_);
	vec_negate(&force_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	/* quadrupole - quadrupole */
	efp_quadrupole_quadrupole_grad(pt_i->quadrupole, pt_j->quadrupole,
				       &dr, &force_, &torque_i_, &torque_j_);
	vec_negate(&torque_j_);
	add_3(&force, &force_, &torque_i, &torque_i_, &torque_j, &torque_j_);

	vec_scale(&force, swf->swf);
	vec_scale(&torque_i, swf->swf);
	vec_scale(&torque_j, swf->swf);

	efp_add_force(efp->grad + fr_i_idx, CVEC(fr_i->x), CVEC(pt_i->x), &force, &torque_i);
	efp_sub_force(efp->grad + fr_j_idx, CVEC(fr_j->x), CVEC(pt_j->x), &force, &torque_j);
	efp_add_stress(&swf->dr, &force, &efp->stress);
}

double
efp_frag_frag_elec(struct efp *efp, size_t fr_i_idx, size_t fr_j_idx)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	double energy = 0.0;
	struct swf swf = efp_make_swf(efp, fr_i, fr_j);

	/* nuclei - nuclei */
	for (size_t ii = 0; ii < fr_i->n_atoms; ii++) {
		for (size_t jj = 0; jj < fr_j->n_atoms; jj++) {
			struct efp_atom *at_i = fr_i->atoms + ii;
			struct efp_atom *at_j = fr_j->atoms + jj;

			vec_t dr = {
				at_j->x - at_i->x - swf.cell.x,
				at_j->y - at_i->y - swf.cell.y,
				at_j->z - at_i->z - swf.cell.z
			};

			energy += efp_charge_charge_energy(at_i->znuc, at_j->znuc, &dr);

			if (efp->do_gradient) {
				vec_t force, add_i, add_j;

				efp_charge_charge_grad(at_i->znuc, at_j->znuc,
						       &dr, &force, &add_i, &add_j);
				vec_scale(&force, swf.swf);
				efp_add_force(efp->grad + fr_i_idx, CVEC(fr_i->x),
						CVEC(at_i->x), &force, NULL);
				efp_sub_force(efp->grad + fr_j_idx, CVEC(fr_j->x),
						CVEC(at_j->x), &force, NULL);
				efp_add_stress(&swf.dr, &force, &efp->stress);
			}
		}
	}

	/* nuclei - mult points */
	for (size_t ii = 0; ii < fr_i->n_atoms; ii++) {
		for (size_t jj = 0; jj < fr_j->n_multipole_pts; jj++) {
			energy += atom_mult_energy(efp, fr_i, fr_j, ii, jj, &swf);

			if (efp->do_gradient)
				atom_mult_grad(efp, fr_i_idx, fr_j_idx, ii, jj, &swf);
		}
	}

	/* mult points - nuclei */
	for (size_t jj = 0; jj < fr_j->n_atoms; jj++) {
		for (size_t ii = 0; ii < fr_i->n_multipole_pts; ii++) {
			struct swf swf2 = swf;

			vec_negate(&swf2.cell);
			vec_negate(&swf2.dr);
			vec_negate(&swf2.dswf);

			energy += atom_mult_energy(efp, fr_j, fr_i, jj, ii, &swf2);

			if (efp->do_gradient)
				atom_mult_grad(efp, fr_j_idx, fr_i_idx, jj, ii, &swf2);
		}
	}

	/* mult points - mult points */
	for (size_t ii = 0; ii < fr_i->n_multipole_pts; ii++) {
		for (size_t jj = 0; jj < fr_j->n_multipole_pts; jj++) {
			energy += mult_mult_energy(efp, fr_i_idx, fr_j_idx, ii, jj, &swf);

			if (efp->do_gradient)
				mult_mult_grad(efp, fr_i_idx, fr_j_idx, ii, jj, &swf);
		}
	}

	vec_t force = {
		swf.dswf.x * energy,
		swf.dswf.y * energy,
		swf.dswf.z * energy
	};

	six_atomic_add_xyz(efp->grad + fr_i_idx, &force);
	six_atomic_sub_xyz(efp->grad + fr_j_idx, &force);
	efp_add_stress(&swf.dr, &force, &efp->stress);

	return energy * swf.swf;
}

static void
rotate_quadrupole(const mat_t *rotmat, const double *in, double *out)
{
	double full_in[9], full_out[9];

	for (size_t a = 0; a < 3; a++)
		for (size_t b = 0; b < 3; b++)
			full_in[a * 3 + b] = in[quad_idx(a, b)];

	efp_rotate_t2(rotmat, full_in, full_out);

	for (size_t a = 0; a < 3; a++)
		for (size_t b = 0; b < 3; b++)
			out[quad_idx(a, b)] = full_out[a * 3 + b];
}

static void
rotate_octupole(const mat_t *rotmat, const double *in, double *out)
{
	double full_in[27], full_out[27];

	for (size_t a = 0; a < 3; a++)
		for (size_t b = 0; b < 3; b++)
			for (size_t c = 0; c < 3; c++) {
				size_t idx = 9 * a + 3 * b + c;
				full_in[idx] = in[oct_idx(a, b, c)];
			}

	efp_rotate_t3(rotmat, full_in, full_out);

	for (size_t a = 0; a < 3; a++)
		for (size_t b = 0; b < 3; b++)
			for (size_t c = 0; c < 3; c++) {
				size_t idx = 9 * a + 3 * b + c;
				out[oct_idx(a, b, c)] = full_out[idx];
			}
}

void
efp_update_elec(struct frag *frag)
{
	for (size_t i = 0; i < frag->n_multipole_pts; i++) {
		const struct multipole_pt *in = frag->lib->multipole_pts + i;
		struct multipole_pt *out = frag->multipole_pts + i;

		/* move point position */
		efp_move_pt(CVEC(frag->x), &frag->rotmat, CVEC(in->x), VEC(out->x));

		/* rotate dipole */
		out->dipole = mat_vec(&frag->rotmat, &in->dipole);

		/* rotate quadrupole */
		rotate_quadrupole(&frag->rotmat, in->quadrupole, out->quadrupole);

		/* correction for Buckingham quadrupoles */
		double *quad = out->quadrupole;

		double qtr = quad[quad_idx(0, 0)] +
			     quad[quad_idx(1, 1)] +
			     quad[quad_idx(2, 2)];

		quad[0] = 1.5 * quad[0] - 0.5 * qtr;
		quad[1] = 1.5 * quad[1] - 0.5 * qtr;
		quad[2] = 1.5 * quad[2] - 0.5 * qtr;
		quad[3] = 1.5 * quad[3];
		quad[4] = 1.5 * quad[4];
		quad[5] = 1.5 * quad[5];

		/* rotate octupole */
		rotate_octupole(&frag->rotmat, in->octupole, out->octupole);

		/* correction for Buckingham octupoles */
		double *oct = out->octupole;

		double otrx = oct[oct_idx(0, 0, 0)] +
			      oct[oct_idx(0, 1, 1)] +
			      oct[oct_idx(0, 2, 2)];
		double otry = oct[oct_idx(0, 0, 1)] +
			      oct[oct_idx(1, 1, 1)] +
			      oct[oct_idx(1, 2, 2)];
		double otrz = oct[oct_idx(0, 0, 2)] +
			      oct[oct_idx(1, 1, 2)] +
			      oct[oct_idx(2, 2, 2)];

		oct[0] = 2.5 * oct[0] - 1.5 * otrx;
		oct[1] = 2.5 * oct[1] - 1.5 * otry;
		oct[2] = 2.5 * oct[2] - 1.5 * otrz;
		oct[3] = 2.5 * oct[3] - 0.5 * otry;
		oct[4] = 2.5 * oct[4] - 0.5 * otrz;
		oct[5] = 2.5 * oct[5] - 0.5 * otrx;
		oct[6] = 2.5 * oct[6] - 0.5 * otrz;
		oct[7] = 2.5 * oct[7] - 0.5 * otrx;
		oct[8] = 2.5 * oct[8] - 0.5 * otry;
		oct[9] = 2.5 * oct[9];
	}
}

static double
compute_ai_elec_frag(struct efp *efp, size_t frag_idx)
{
	double energy = 0.0;
	struct frag *fr_i = efp->frags + frag_idx;

	for (size_t i = 0; i < fr_i->n_atoms; i++) {
		for (size_t j = 0; j < efp->n_ptc; j++) {
			struct efp_atom *at_i = fr_i->atoms + i;
			vec_t dr = vec_sub(CVEC(at_i->x), efp->ptc_xyz + j);

			energy += efp_charge_charge_energy(at_i->znuc, efp->ptc[j], &dr);
		}
	}
	for (size_t i = 0; i < fr_i->n_multipole_pts; i++) {
		for (size_t j = 0; j < efp->n_ptc; j++) {
			struct multipole_pt *pt_i = fr_i->multipole_pts + i;
			vec_t dr = vec_sub(CVEC(pt_i->x), efp->ptc_xyz + j);

			/* charge - monopole */
			energy += efp_charge_charge_energy(efp->ptc[j], pt_i->monopole, &dr);

			/* charge - dipole */
			energy += efp_charge_dipole_energy(efp->ptc[j], &pt_i->dipole, &dr);

			/* charge - quadrupole */
			energy += efp_charge_quadrupole_energy(efp->ptc[j],
					pt_i->quadrupole, &dr);

			/* charge - octupole */
			energy += efp_charge_octupole_energy(efp->ptc[j],
					pt_i->octupole, &dr);
		}
	}
	return energy;
}

static void
compute_ai_elec_frag_grad(struct efp *efp, size_t frag_idx)
{
	struct frag *fr_j = efp->frags + frag_idx;

	for (size_t i = 0; i < efp->n_ptc; i++) {
		/* ab initio atom - fragment atoms */
		for (size_t k = 0; k < fr_j->n_atoms; k++) {
			struct efp_atom *at_j = fr_j->atoms + k;

			vec_t force, add_i, add_j;
			vec_t dr = vec_sub(CVEC(at_j->x), efp->ptc_xyz + i);

			efp_charge_charge_grad(efp->ptc[i], at_j->znuc, &dr,
					       &force, &add_i, &add_j);

			vec_atomic_add(efp->ptc_grad + i, &force);
			efp_sub_force(efp->grad + frag_idx, CVEC(fr_j->x), CVEC(at_j->x),
					&force, &add_j);
		}

		/* ab initio atom - fragment multipoles */
		for (size_t k = 0; k < fr_j->n_multipole_pts; k++) {
			struct multipole_pt *pt_j = fr_j->multipole_pts + k;

			vec_t force_, add_i_, add_j_;
			vec_t force = vec_zero, add_i = vec_zero, add_j = vec_zero;
			vec_t dr = vec_sub(CVEC(pt_j->x), efp->ptc_xyz + i);

			/* monopole */
			efp_charge_charge_grad(efp->ptc[i], pt_j->monopole, &dr,
					       &force_, &add_i_, &add_j_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			/* dipole */
			efp_charge_dipole_grad(efp->ptc[i], &pt_j->dipole, &dr,
					       &force_, &add_i_, &add_j_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			/* quadrupole */
			efp_charge_quadrupole_grad(efp->ptc[i], pt_j->quadrupole, &dr,
						   &force_, &add_i_, &add_j_);
			vec_negate(&add_j_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			/* octupole */
			efp_charge_octupole_grad(efp->ptc[i], pt_j->octupole, &dr,
						 &force_, &add_i_, &add_j_);
			add_3(&force, &force_, &add_i, &add_i_, &add_j, &add_j_);

			vec_atomic_add(efp->ptc_grad + i, &force);
			efp_sub_force(efp->grad + frag_idx, CVEC(fr_j->x), CVEC(pt_j->x),
					&force, &add_j);
		}
	}
}

static void
compute_ai_elec_range(struct efp *efp, size_t from, size_t to, void *data)
{
	double energy = 0.0;

	(void)data;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:energy)
#endif
	for (size_t i = from; i < to; i++) {
		energy += compute_ai_elec_frag(efp, i);

		if (efp->do_gradient)
			compute_ai_elec_frag_grad(efp, i);
	}

	efp->energy.electrostatic_point_charges += energy;
}

enum efp_result
efp_compute_ai_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_AI_ELEC))
		return (EFP_RESULT_SUCCESS);

	efp_balance_work(efp, compute_ai_elec_range, NULL);
	efp_allreduce(&efp->energy.electrostatic_point_charges, 1);

	return (EFP_RESULT_SUCCESS);
}
