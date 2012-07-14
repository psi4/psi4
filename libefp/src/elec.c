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

#include "efp_private.h"
#include "elec.h"

static double
get_screen_damping(double r_ij, double pi, double pj)
{
	if (pj == HUGE_VAL) {   /* j is nucleus */
		return 1.0 - exp(-pi * r_ij);
	}
	else if (fabs((pi - pj) * r_ij) < 1.0e-5) {
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
	else if (fabs((pi - pj) * r_ij) < 1.0e-5) {
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
charge_charge_energy(double q1, double q2, double r)
{
	return q1 * q2 / r;
}

static double
charge_dipole_energy(double q1, const vec_t *d2, const vec_t *dr, double r)
{
	double r3 = r * r * r;

	return -q1 / r3 * vec_dot(d2, dr);
}

static double
charge_quadrupole_energy(double q1, const double *quad2, const vec_t *dr,
			 double r)
{
	double r2 = r * r;
	double r5 = r2 * r2 * r;

	return q1 / r5 * quadrupole_sum(quad2, dr);
}

static double
charge_octupole_energy(double q1, const double *oct2, const vec_t *dr,
		       double r)
{
	double r2 = r * r;
	double r7 = r2 * r2 * r2 * r;

	return -q1 / r7 * octupole_sum(oct2, dr);
}

static double
dipole_dipole_energy(const vec_t *d1, const vec_t *d2, const vec_t *dr,
		     double r)
{
	double r2 = r * r;
	double r3 = r2 * r;
	double r5 = r3 * r2;

	double d1dr = vec_dot(d1, dr);
	double d2dr = vec_dot(d2, dr);

	return vec_dot(d1, d2) / r3 - 3.0 * d1dr * d2dr / r5;
}

static double
dipole_quadrupole_energy(const vec_t *d1, const double *quad2, const vec_t *dr,
			 double r)
{
	double r2 = r * r;
	double r5 = r2 * r2 * r;
	double r7 = r5 * r2;

	double d1dr = vec_dot(d1, dr);
	double q2dr = quadrupole_sum(quad2, dr);
	double d1q2dr = 0.0;

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++) {
			int idx = quad_idx(a, b);
			d1q2dr += quad2[idx] * vec_el(d1, a) * vec_el(dr, b);
		}

	return 5.0 / r7 * q2dr * d1dr - 2.0 / r5 * d1q2dr;
}

static double
quadrupole_quadrupole_energy(const double *quad1, const double *quad2,
			     const vec_t *dr, double r)
{
	double r2 = r * r;
	double r5 = r2 * r2 * r;
	double r7 = r5 * r2;
	double r9 = r7 * r2;

	double q1dr = quadrupole_sum(quad1, dr);
	double q2dr = quadrupole_sum(quad2, dr);

	double q1q2 = 0.0;
	double q1q2dr = 0.0;

	for (int a = 0; a < 3; a++) {
		double t1 = 0.0;
		double t2 = 0.0;

		for (int b = 0; b < 3; b++) {
			int idx = quad_idx(a, b);

			t1 += quad1[idx] * vec_el(dr, b);
			t2 += quad2[idx] * vec_el(dr, b);

			q1q2 += quad1[idx] * quad2[idx];
		}

		q1q2dr += t1 * t2;
	}

	return (2.0 / r5 * q1q2 - 20.0 / r7 * q1q2dr +
						35.0 / r9 * q1dr * q2dr) / 3.0;
}

static double
atom_mult_energy(struct efp *efp, double charge, const vec_t *pos,
		 int frag_idx, int pt_idx)
{
	struct frag *fr_i = efp->frags + frag_idx;
	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_idx;

	vec_t dr = vec_sub(VEC(pt_i->x), pos);
	double r = vec_len(&dr);

	double energy = 0.0, ccdamp = 1.0;

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double sp = fr_i->screen_params[pt_idx];
		ccdamp = get_screen_damping(r, sp, HUGE_VAL);
	}

	/* charge - monopole */
	energy += ccdamp * charge_charge_energy(
				charge, pt_i->monopole, r);

	/* charge - dipole */
	energy += charge_dipole_energy(
				charge, &pt_i->dipole, &dr, r);

	/* charge - quadrupole */
	energy += charge_quadrupole_energy(
				charge, pt_i->quadrupole, &dr, r);

	/* charge - octupole */
	energy += charge_octupole_energy(
				charge, pt_i->octupole, &dr, r);

	return energy;
}

static void
atom_mult_grad(struct efp *efp, int fr_i_idx, int fr_j_idx,
	       int atom_i_idx, int pt_j_idx)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	struct efp_atom *at_i = fr_i->atoms + atom_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = vec_sub(VEC(pt_j->x), VEC(at_i->x));
	vec_t force, torque_i, torque_j;

	/* charge - charge */
	efp_charge_charge_grad(at_i->znuc, pt_j->monopole, &dr,
			       &force, &torque_i, &torque_j);

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double r = vec_len(&dr);
		double sp = fr_j->screen_params[pt_j_idx];
		double gdamp = get_screen_damping_grad(r, sp, HUGE_VAL);

		force.x *= gdamp;
		force.y *= gdamp;
		force.z *= gdamp;
	}

	add_force_torque_2(fr_i, fr_j, VEC(at_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* charge - dipole */
	efp_charge_dipole_grad(at_i->znuc, &pt_j->dipole, &dr,
			       &force, &torque_i, &torque_j);
	add_force_torque_2(fr_i, fr_j, VEC(at_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* charge - quadrupole */
	efp_charge_quadrupole_grad(at_i->znuc, pt_j->quadrupole, &dr,
				   &force, &torque_i, &torque_j);

	vec_negate(&torque_j);

	add_force_torque_2(fr_i, fr_j, VEC(at_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* charge - octupole */
	efp_charge_octupole_grad(at_i->znuc, pt_j->octupole, &dr,
				 &force, &torque_i, &torque_j);
	add_force_torque_2(fr_i, fr_j, VEC(at_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);
}

static double
mult_mult_energy(struct efp *efp, int fr_i_idx, int fr_j_idx,
		 int pt_i_idx, int pt_j_idx)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = vec_sub(VEC(pt_j->x), VEC(pt_i->x));
	double r = vec_len(&dr);

	double energy = 0.0, ccdamp = 1.0;

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double screen_i = fr_i->screen_params[pt_i_idx];
		double screen_j = fr_j->screen_params[pt_j_idx];
		ccdamp = get_screen_damping(r, screen_i, screen_j);
	}

	/* monopole - monopole */
	energy += ccdamp * charge_charge_energy(
				pt_i->monopole, pt_j->monopole, r);

	/* monopole - dipole */
	energy += charge_dipole_energy(
				pt_i->monopole, &pt_j->dipole, &dr, r);

	/* dipole - monopole */
	energy -= charge_dipole_energy(
				pt_j->monopole, &pt_i->dipole, &dr, r);

	/* monopole - quadrupole */
	energy += charge_quadrupole_energy(
				pt_i->monopole, pt_j->quadrupole, &dr, r);

	/* quadrupole - monopole */
	energy += charge_quadrupole_energy(
				pt_j->monopole, pt_i->quadrupole, &dr, r);

	/* monopole - octupole */
	energy += charge_octupole_energy(
				pt_i->monopole, pt_j->octupole, &dr, r);

	/* octupole - monopole */
	energy -= charge_octupole_energy(
				pt_j->monopole, pt_i->octupole, &dr, r);

	/* dipole - dipole */
	energy += dipole_dipole_energy(
				&pt_i->dipole, &pt_j->dipole, &dr, r);

	/* dipole - quadrupole */
	energy += dipole_quadrupole_energy(
				&pt_i->dipole, pt_j->quadrupole, &dr, r);

	/* quadrupole - dipole */
	energy -= dipole_quadrupole_energy(
				&pt_j->dipole, pt_i->quadrupole, &dr, r);

	/* quadrupole - quadrupole */
	energy += quadrupole_quadrupole_energy(
				pt_i->quadrupole, pt_j->quadrupole, &dr, r);

	return energy;
}

static void
mult_mult_grad(struct efp *efp, int fr_i_idx, int fr_j_idx,
	       int pt_i_idx, int pt_j_idx)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = vec_sub(VEC(pt_j->x), VEC(pt_i->x));
	vec_t force, torque_i, torque_j;

	/* monopole - monopole */
	efp_charge_charge_grad(pt_i->monopole, pt_j->monopole, &dr,
			       &force, &torque_i, &torque_j);

	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double r = vec_len(&dr);

		double screen_i = fr_i->screen_params[pt_i_idx];
		double screen_j = fr_j->screen_params[pt_j_idx];

		double gdamp = get_screen_damping_grad(r, screen_i, screen_j);

		force.x *= gdamp;
		force.y *= gdamp;
		force.z *= gdamp;
	}

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* monopole - dipole */
	efp_charge_dipole_grad(pt_i->monopole, &pt_j->dipole, &dr,
			       &force, &torque_i, &torque_j);
	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* dipole - monopole */
	efp_charge_dipole_grad(pt_j->monopole, &pt_i->dipole, &dr,
			       &force, &torque_j, &torque_i);

	vec_negate(&force);

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* monopole - quadrupole */
	efp_charge_quadrupole_grad(pt_i->monopole, pt_j->quadrupole, &dr,
				   &force, &torque_i, &torque_j);

	vec_negate(&torque_j);

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* quadrupole - monopole */
	efp_charge_quadrupole_grad(pt_j->monopole, pt_i->quadrupole, &dr,
				   &force, &torque_j, &torque_i);
	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* monopole - octupole */
	efp_charge_octupole_grad(pt_i->monopole, pt_j->octupole, &dr,
				 &force, &torque_i, &torque_j);
	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* octupole - monopole */
	efp_charge_octupole_grad(pt_j->monopole, pt_i->octupole, &dr,
				 &force, &torque_j, &torque_i);

	vec_negate(&force);

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* dipole - dipole */
	efp_dipole_dipole_grad(&pt_i->dipole, &pt_j->dipole, &dr,
			       &force, &torque_i, &torque_j);

	vec_negate(&torque_j);

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* dipole - quadrupole */
	efp_dipole_quadrupole_grad(&pt_i->dipole, pt_j->quadrupole, &dr,
				   &force, &torque_i, &torque_j);
	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* quadrupole - dipole */
	efp_dipole_quadrupole_grad(&pt_j->dipole, pt_i->quadrupole, &dr,
				   &force, &torque_j, &torque_i);

	vec_negate(&force);

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);

	/* quadrupole - quadrupole */
	efp_quadrupole_quadrupole_grad(pt_i->quadrupole, pt_j->quadrupole,
				       &dr, &force, &torque_i, &torque_j);

	vec_negate(&torque_j);

	add_force_torque_2(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
			   &force, &torque_i, &torque_j);
}

static double
frag_frag_elec(struct efp *efp, int fr_i_idx, int fr_j_idx)
{
	double energy = 0.0;

	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	/* nuclei - nuclei */
	for (int ii = 0; ii < fr_i->n_atoms; ii++) {
		for (int jj = 0; jj < fr_j->n_atoms; jj++) {
			struct efp_atom *at_i = fr_i->atoms + ii;
			struct efp_atom *at_j = fr_j->atoms + jj;

			vec_t dr = vec_sub(VEC(at_j->x), VEC(at_i->x));
			double r = vec_len(&dr);

			energy += charge_charge_energy(at_i->znuc,
						       at_j->znuc, r);

			if (efp->do_gradient) {
				vec_t force, add_i, add_j;

				efp_charge_charge_grad(at_i->znuc, at_j->znuc,
						       &dr, &force, &add_i,
						       &add_j);
				add_force_torque_2(fr_i, fr_j,
						   VEC(at_i->x), VEC(at_j->x),
						   &force, &add_i, &add_j);
			}
		}
	}

	/* nuclei - mult points */
	for (int ii = 0; ii < fr_i->n_atoms; ii++) {
		for (int jj = 0; jj < fr_j->n_multipole_pts; jj++) {
			struct efp_atom *at_i = fr_i->atoms + ii;

			energy += atom_mult_energy(efp, at_i->znuc,
						   VEC(at_i->x), fr_j_idx, jj);

			if (efp->do_gradient)
				atom_mult_grad(efp, fr_i_idx, fr_j_idx, ii, jj);
		}
	}

	/* mult points - nuclei */
	for (int jj = 0; jj < fr_j->n_atoms; jj++) {
		for (int ii = 0; ii < fr_i->n_multipole_pts; ii++) {
			struct efp_atom *at_j = fr_j->atoms + jj;

			energy += atom_mult_energy(efp, at_j->znuc,
						   VEC(at_j->x), fr_i_idx, ii);

			if (efp->do_gradient)
				atom_mult_grad(efp, fr_j_idx, fr_i_idx, jj, ii);
		}
	}

	/* mult points - mult points */
	for (int ii = 0; ii < fr_i->n_multipole_pts; ii++)
		for (int jj = 0; jj < fr_j->n_multipole_pts; jj++) {
			energy += mult_mult_energy(efp, fr_i_idx, fr_j_idx,
						   ii, jj);

			if (efp->do_gradient)
				mult_mult_grad(efp, fr_i_idx, fr_j_idx, ii, jj);
		}

	return energy;
}

enum efp_result
efp_compute_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_ELEC))
		return EFP_RESULT_SUCCESS;

	double energy = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:energy)
	for (int i = 0; i < efp->n_frag; i++)
		for (int j = i + 1; j < efp->n_frag; j++)
			energy += frag_frag_elec(efp, i, j);

	efp->energy.electrostatic = energy;
	return EFP_RESULT_SUCCESS;
}

static void
rotate_quadrupole(const mat_t *rotmat, const double *in, double *out)
{
	double full_in[9], full_out[9];

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			full_in[a * 3 + b] = in[quad_idx(a, b)];

	rotate_t2(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			out[quad_idx(a, b)] = full_out[a * 3 + b];
}

static void
rotate_octupole(const mat_t *rotmat, const double *in, double *out)
{
	double full_in[27], full_out[27];

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int idx = 9 * a + 3 * b + c;
				full_in[idx] = in[oct_idx(a, b, c)];
			}

	rotate_t3(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int idx = 9 * a + 3 * b + c;
				out[oct_idx(a, b, c)] = full_out[idx];
			}
}

void
efp_update_elec(struct frag *frag)
{
	const mat_t *rotmat = &frag->rotmat;

	for (int i = 0; i < frag->n_multipole_pts; i++) {
		const struct multipole_pt *pt_in =
					frag->lib->multipole_pts + i;
		struct multipole_pt *pt_out =
					frag->multipole_pts + i;

		/* move point position */
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			VEC(pt_in->x), VEC(pt_out->x));

		/* rotate dipole */
		mat_vec(rotmat, &pt_in->dipole, &pt_out->dipole);

		/* rotate quadrupole */
		rotate_quadrupole(rotmat, pt_in->quadrupole,
				  pt_out->quadrupole);

		/* correction for Buckingham quadrupoles */
		double *quad = pt_out->quadrupole;

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
		rotate_octupole(rotmat, pt_in->octupole, pt_out->octupole);

		/* correction for Buckingham octupoles */
		double *oct = pt_out->octupole;

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
compute_ai_elec_frag(struct efp *efp, int frag_idx)
{
	double energy = 0.0;
	struct frag *fr_i = efp->frags + frag_idx;

	for (int i = 0; i < fr_i->n_atoms; i++) {
		for (int j = 0; j < efp->qm.n_atoms; j++) {
			struct efp_atom *at_i = fr_i->atoms + i;
			const double *xyz_j = efp->qm.xyz + 3 * j;
			double znuc_j = efp->qm.znuc[j];

			double r = vec_dist(VEC(at_i->x), (const vec_t *)xyz_j);
			energy += at_i->znuc * znuc_j / r;
		}
	}
	for (int i = 0; i < fr_i->n_multipole_pts; i++) {
		for (int j = 0; j < efp->qm.n_atoms; j++) {
			const double *xyz = efp->qm.xyz + 3 * j;
			double znuc = efp->qm.znuc[j];

			energy += atom_mult_energy(efp, znuc,
					(const vec_t *)xyz, frag_idx, i);
		}
	}
	return energy;
}

static void
compute_ai_elec_frag_grad(struct efp *efp, int frag_idx)
{
	/* XXX */
	assert(0);
}

enum efp_result
efp_compute_ai_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_AI_ELEC))
		return EFP_RESULT_SUCCESS;

	double energy = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:energy)
	for (int i = 0; i < efp->n_frag; i++) {
		energy += compute_ai_elec_frag(efp, i);

		if (efp->do_gradient)
			compute_ai_elec_frag_grad(efp, i);
	}

	efp->energy.ai_electrostatic = energy;
	return EFP_RESULT_SUCCESS;
}
