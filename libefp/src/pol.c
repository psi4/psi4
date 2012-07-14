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

static void
add_multipole_field(struct polarizable_pt *pt,
		    const struct multipole_pt *mult_pt)
{
	vec_t dr = vec_sub(VEC(pt->x), VEC(mult_pt->x));

	double t1, t2;

	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;
	double r7 = r5 * r * r;

	/* charge */
	#pragma omp atomic
	pt->elec_field.x += mult_pt->monopole * dr.x / r3;
	#pragma omp atomic
	pt->elec_field.y += mult_pt->monopole * dr.y / r3;
	#pragma omp atomic
	pt->elec_field.z += mult_pt->monopole * dr.z / r3;

	/* dipole */
	t1 = vec_dot(&mult_pt->dipole, &dr);

	#pragma omp atomic
	pt->elec_field.x -= mult_pt->dipole.x / r3 - 3.0 / r5 * t1 * dr.x;
	#pragma omp atomic
	pt->elec_field.y -= mult_pt->dipole.y / r3 - 3.0 / r5 * t1 * dr.y;
	#pragma omp atomic
	pt->elec_field.z -= mult_pt->dipole.z / r3 - 3.0 / r5 * t1 * dr.z;

	/* quadrupole */
	t1 = quadrupole_sum(mult_pt->quadrupole, &dr);

	t2 = mult_pt->quadrupole[quad_idx(0, 0)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 0)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 0)] * dr.z;
	#pragma omp atomic
	pt->elec_field.x -= 2.0 / r5 * t2 - 5.0 / r7 * t1 * dr.x;

	t2 = mult_pt->quadrupole[quad_idx(0, 1)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 1)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 1)] * dr.z;
	#pragma omp atomic
	pt->elec_field.y -= 2.0 / r5 * t2 - 5.0 / r7 * t1 * dr.y;

	t2 = mult_pt->quadrupole[quad_idx(0, 2)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 2)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 2)] * dr.z;
	#pragma omp atomic
	pt->elec_field.z -= 2.0 / r5 * t2 - 5.0 / r7 * t1 * dr.z;

	/* octupole-polarizability interactions are ignored */
}

static void
compute_elec_field_pt(struct efp *efp, int frag_idx, int pt_idx)
{
	struct frag *frag = efp->frags + frag_idx;
	struct polarizable_pt *pt = frag->polarizable_pts + pt_idx;

	vec_zero(&pt->elec_field);

	for (int i = 0; i < efp->n_frag; i++) {
		if (i == frag_idx)
			continue;

		/* field due to nuclei */
		for (int j = 0; j < efp->frags[i].n_atoms; j++) {
			struct efp_atom *at = efp->frags[i].atoms + j;

			vec_t dr = vec_sub(VEC(pt->x), VEC(at->x));

			double r = vec_len(&dr);
			double r3 = r * r * r;

			#pragma omp atomic
			pt->elec_field.x += at->znuc * dr.x / r3;
			#pragma omp atomic
			pt->elec_field.y += at->znuc * dr.y / r3;
			#pragma omp atomic
			pt->elec_field.z += at->znuc * dr.z / r3;
		}

		/* field due to multipoles */
		for (int j = 0; j < efp->frags[i].n_multipole_pts; j++) {
			struct multipole_pt *mult_pt =
					efp->frags[i].multipole_pts + j;
			add_multipole_field(pt, mult_pt);
		}
	}

	if (efp->opts.terms & EFP_TERM_AI_POL) {
		/* field due to QM nuclei */
		for (int i = 0; i < efp->qm.n_atoms; i++) {
			const double *xyz = efp->qm.xyz + 3 * i;
			double znuc = efp->qm.znuc[i];

			vec_t dr = vec_sub(VEC(pt->x), (const vec_t *)xyz);

			double r = vec_len(&dr);
			double r3 = r * r * r;

			#pragma omp atomic
			pt->elec_field.x += znuc * dr.x / r3;
			#pragma omp atomic
			pt->elec_field.y += znuc * dr.y / r3;
			#pragma omp atomic
			pt->elec_field.z += znuc * dr.z / r3;
		}
	}
}

static void
add_electron_density_field(struct efp *efp)
{
	int n_pt = 0;
	for (int i = 0; i < efp->n_frag; i++)
		n_pt += efp->frags[i].n_polarizable_pts;

	double xyz[3 * n_pt], field[3 * n_pt], *ptr;

	ptr = xyz;
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			*ptr++ = pt->x;
			*ptr++ = pt->y;
			*ptr++ = pt->z;
		}
	}

	efp->callbacks.get_electron_density_field(n_pt, xyz, field,
			efp->callbacks.get_electron_density_field_user_data);

	ptr = field;
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			pt->elec_field.x += *ptr++;
			pt->elec_field.y += *ptr++;
			pt->elec_field.z += *ptr++;
		}
	}
}

static void
compute_elec_field(struct efp *efp)
{
	#pragma omp parallel for schedule(dynamic, 4)
	for (int i = 0; i < efp->n_frag; i++)
		for (int j = 0; j < efp->frags[i].n_polarizable_pts; j++)
			compute_elec_field_pt(efp, i, j);

	if (efp->opts.terms & EFP_TERM_AI_POL)
		add_electron_density_field(efp);
}

static void
get_induced_dipole_field(struct efp *efp, int frag_idx,
			 struct polarizable_pt *pt,
			 vec_t *field, vec_t *field_conj)
{
	vec_zero(field);
	vec_zero(field_conj);

	for (int j = 0; j < efp->n_frag; j++) {
		if (j == frag_idx)
			continue;

		struct frag *fr_j = efp->frags + j;

		for (int jj = 0; jj < fr_j->n_polarizable_pts; jj++) {
			struct polarizable_pt *pt_j =
				fr_j->polarizable_pts + jj;

			vec_t dr = vec_sub(VEC(pt->x), VEC(pt_j->x));

			double ir = 1.0 / vec_len(&dr);
			double ir3 = ir * ir * ir;
			double ir5 = ir3 * ir * ir;

			double t1 = vec_dot(&pt_j->induced_dipole, &dr);

			field->x -= pt_j->induced_dipole.x * ir3 -
						3.0 * t1 * dr.x * ir5;
			field->y -= pt_j->induced_dipole.y * ir3 -
						3.0 * t1 * dr.y * ir5;
			field->z -= pt_j->induced_dipole.z * ir3 -
						3.0 * t1 * dr.z * ir5;

			double t2 = vec_dot(&pt_j->induced_dipole_conj, &dr);

			field_conj->x -= pt_j->induced_dipole_conj.x * ir3 -
						3.0 * t2 * dr.x * ir5;
			field_conj->y -= pt_j->induced_dipole_conj.y * ir3 -
						3.0 * t2 * dr.y * ir5;
			field_conj->z -= pt_j->induced_dipole_conj.z * ir3 -
						3.0 * t2 * dr.z * ir5;
		}
	}
}

static double
pol_scf_iter(struct efp *efp)
{
	/* compute new induced dipoles on polarizable points */
	#pragma omp parallel for schedule(dynamic, 4)
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			/* electric field from other induced dipoles */
			vec_t field, field_conj;
			get_induced_dipole_field(efp, i, pt, &field,
					&field_conj);

			/* add field that doesn't change during scf */
			field.x += pt->elec_field.x;
			field.y += pt->elec_field.y;
			field.z += pt->elec_field.z;

			field_conj.x += pt->elec_field.x;
			field_conj.y += pt->elec_field.y;
			field_conj.z += pt->elec_field.z;

			mat_vec(&pt->tensor, &field, &pt->induced_dipole_new);
			mat_trans_vec(&pt->tensor, &field_conj,
					&pt->induced_dipole_conj_new);
		}
	}

	int n_pt = 0;
	double conv = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:n_pt,conv)
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		n_pt += frag->n_polarizable_pts;

		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			conv += vec_dist_2(&pt->induced_dipole_new,
					   &pt->induced_dipole);

			pt->induced_dipole = pt->induced_dipole_new;
			pt->induced_dipole_conj = pt->induced_dipole_conj_new;
		}
	}

	return conv / n_pt;
}

double
efp_compute_pol_energy(struct efp *efp)
{
	static const double scf_conv_epsilon = 1.0e-16;

	compute_elec_field(efp);

	/* set initial approximation - all induced dipoles are zero */

	#pragma omp parallel for schedule(dynamic, 4)
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			vec_zero(&pt->induced_dipole);
			vec_zero(&pt->induced_dipole_conj);
		}
	}

	/* compute induced dipoles self consistently */
	while (pol_scf_iter(efp) > scf_conv_epsilon);

	double energy = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:energy)
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			energy += -0.5 * vec_dot(&pt->induced_dipole,
						 &pt->elec_field);
		}
	}
	return energy;
}

static void
compute_grad_point(struct efp *efp, int frag_idx, int pt_idx)
{
	struct frag *fr_i = efp->frags + frag_idx;
	struct polarizable_pt *pt_i = fr_i->polarizable_pts + pt_idx;

	vec_t dipole_i = {
		0.5 * (pt_i->induced_dipole.x + pt_i->induced_dipole_conj.x),
		0.5 * (pt_i->induced_dipole.y + pt_i->induced_dipole_conj.y),
		0.5 * (pt_i->induced_dipole.z + pt_i->induced_dipole_conj.z)
	};

	for (int j = 0; j < efp->n_frag; j++) {
		if (j == frag_idx)
			continue;

		struct frag *fr_j = efp->frags + j;

		/* induced dipole - nuclei */
		for (int k = 0; k < fr_j->n_atoms; k++) {
			struct efp_atom *at_j = fr_j->atoms + k;

			vec_t dr = vec_sub(VEC(at_j->x), VEC(pt_i->x));
			vec_t force, add_i, add_j;

			efp_charge_dipole_grad(at_j->znuc, &dipole_i, &dr,
					       &force, &add_j, &add_i);
			vec_negate(&force);

			add_force_torque_2(fr_i, fr_j,
					   VEC(pt_i->x), VEC(at_j->x),
					   &force, &add_i, &add_j);
		}

		/* induced dipole - multipoles */
		for (int k = 0; k < fr_j->n_multipole_pts; k++) {
			struct multipole_pt *pt_j = fr_j->multipole_pts + k;

			vec_t dr = vec_sub(VEC(pt_j->x), VEC(pt_i->x));
			vec_t force, add_i, add_j;

			/* induced dipole - charge */
			efp_charge_dipole_grad(pt_j->monopole, &dipole_i, &dr,
					       &force, &add_j, &add_i);
			vec_negate(&force);

			add_force_torque_2(fr_i, fr_j,
					   VEC(pt_i->x), VEC(pt_j->x),
					   &force, &add_i, &add_j);

			/* induced dipole - dipole */
			efp_dipole_dipole_grad(&dipole_i, &pt_j->dipole, &dr,
					       &force, &add_i, &add_j);
			vec_negate(&add_j);

			add_force_torque_2(fr_i, fr_j,
					   VEC(pt_i->x), VEC(pt_j->x),
					   &force, &add_i, &add_j);

			/* induced dipole - quadrupole */
			efp_dipole_quadrupole_grad(&dipole_i, pt_j->quadrupole,
						   &dr, &force, &add_i, &add_j);
			add_force_torque_2(fr_i, fr_j,
					   VEC(pt_i->x), VEC(pt_j->x),
					   &force, &add_i, &add_j);

			/* octupole-polarizability interactions are ignored */
		}

		/* induced dipole - induced dipoles */
		for (int k = 0; k < fr_j->n_polarizable_pts; k++) {
			struct polarizable_pt *pt_j =
					fr_j->polarizable_pts + k;

			vec_t dr = vec_sub(VEC(pt_j->x), VEC(pt_i->x));

			vec_t half_dipole_i = {
				0.5 * pt_i->induced_dipole.x,
				0.5 * pt_i->induced_dipole.y,
				0.5 * pt_i->induced_dipole.z
			};

			vec_t force, add_i, add_j;

			efp_dipole_dipole_grad(&half_dipole_i,
					       &pt_j->induced_dipole_conj,
					       &dr, &force, &add_i, &add_j);
			vec_negate(&add_j);

			add_force_torque_2(fr_i, fr_j,
					   VEC(pt_i->x), VEC(pt_j->x),
					   &force, &add_i, &add_j);
		}
	}

	if (efp->opts.terms & EFP_TERM_AI_POL) {
		/* XXX */
		assert(0);
	}
}

static void
compute_grad(struct efp *efp)
{
	#pragma omp parallel for schedule(dynamic, 4)
	for (int i = 0; i < efp->n_frag; i++)
		for (int j = 0; j < efp->frags[i].n_polarizable_pts; j++)
			compute_grad_point(efp, i, j);
}

enum efp_result
efp_compute_pol(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_POL))
		return EFP_RESULT_SUCCESS;

	efp->energy.polarization = efp_compute_pol_energy(efp);

	if (efp->do_gradient)
		compute_grad(efp);

	return EFP_RESULT_SUCCESS;
}

void
efp_update_pol(struct frag *frag)
{
	const mat_t *rotmat = &frag->rotmat;

	for (int i = 0; i < frag->n_polarizable_pts; i++) {
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			VEC(frag->lib->polarizable_pts[i].x),
			VEC(frag->polarizable_pts[i].x));

		const mat_t *in = &frag->lib->polarizable_pts[i].tensor;
		mat_t *out = &frag->polarizable_pts[i].tensor;

		rotate_t2(rotmat, (const double *)in, (double *)out);
	}
}
