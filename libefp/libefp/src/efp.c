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

#include <stdio.h>
#include <stdlib.h>

#include "balance.h"
#include "compat.h"
#include "elec.h"
#include "private.h"
#include "stream.h"

static void
update_fragment(struct frag *frag)
{
	/* update atoms */
	for (size_t i = 0; i < frag->n_atoms; i++)
		efp_move_pt(CVEC(frag->x), &frag->rotmat,
			CVEC(frag->lib->atoms[i].x), VEC(frag->atoms[i].x));

	efp_update_elec(frag);
	efp_update_pol(frag);
	efp_update_disp(frag);
	efp_update_xr(frag);
}

static enum efp_result
set_coord_xyzabc(struct frag *frag, const double *coord)
{
	frag->x = coord[0];
	frag->y = coord[1];
	frag->z = coord[2];

	euler_to_matrix(coord[3], coord[4], coord[5], &frag->rotmat);

	update_fragment(frag);
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
set_coord_points(struct frag *frag, const double *coord)
{
	if (frag->n_atoms < 3) {
		efp_log("fragment must contain at least three atoms");
		return EFP_RESULT_FATAL;
	}

	double ref[9] = {
		frag->lib->atoms[0].x, frag->lib->atoms[0].y, frag->lib->atoms[0].z,
		frag->lib->atoms[1].x, frag->lib->atoms[1].y, frag->lib->atoms[1].z,
		frag->lib->atoms[2].x, frag->lib->atoms[2].y, frag->lib->atoms[2].z
	};

	vec_t p1;
	mat_t rot1, rot2;

	efp_points_to_matrix(coord, &rot1);
	efp_points_to_matrix(ref, &rot2);
	rot2 = mat_transpose(&rot2);
	frag->rotmat = mat_mat(&rot1, &rot2);
	p1 = mat_vec(&frag->rotmat, VEC(frag->lib->atoms[0].x));

	/* center of mass */
	frag->x = coord[0] - p1.x;
	frag->y = coord[1] - p1.y;
	frag->z = coord[2] - p1.z;

	update_fragment(frag);
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
set_coord_rotmat(struct frag *frag, const double *coord)
{
	if (!efp_check_rotation_matrix((const mat_t *)(coord + 3))) {
		efp_log("invalid rotation matrix specified");
		return EFP_RESULT_FATAL;
	}

	frag->x = coord[0];
	frag->y = coord[1];
	frag->z = coord[2];

	memcpy(&frag->rotmat, coord + 3, sizeof(mat_t));

	update_fragment(frag);
	return EFP_RESULT_SUCCESS;
}

static void
free_frag(struct frag *frag)
{
	if (!frag)
		return;

	free(frag->atoms);
	free(frag->multipole_pts);
	free(frag->polarizable_pts);
	free(frag->dynamic_polarizable_pts);
	free(frag->lmo_centroids);
	free(frag->xr_fock_mat);
	free(frag->xr_wf);
	free(frag->xrfit);
	free(frag->screen_params);
	free(frag->ai_screen_params);

	for (size_t i = 0; i < 3; i++)
		free(frag->xr_wf_deriv[i]);

	for (size_t i = 0; i < frag->n_xr_atoms; i++) {
		for (size_t j = 0; j < frag->xr_atoms[i].n_shells; j++)
			free(frag->xr_atoms[i].shells[j].coef);
		free(frag->xr_atoms[i].shells);
	}

	free(frag->xr_atoms);

	/* don't do free(frag) here */
}

static enum efp_result
copy_frag(struct frag *dest, const struct frag *src)
{
	size_t size;
	memcpy(dest, src, sizeof(struct frag));

	if (src->atoms) {
		size = src->n_atoms * sizeof(struct efp_atom);
		dest->atoms = (struct efp_atom *)malloc(size);
		if (!dest->atoms)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->atoms, src->atoms, size);
	}
	if (src->multipole_pts) {
		size = src->n_multipole_pts * sizeof(struct multipole_pt);
		dest->multipole_pts = (struct multipole_pt *)malloc(size);
		if (!dest->multipole_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->multipole_pts, src->multipole_pts, size);
	}
	if (src->screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->screen_params = (double *)malloc(size);
		if (!dest->screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->screen_params, src->screen_params, size);
	}
	if (src->ai_screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->ai_screen_params = (double *)malloc(size);
		if (!dest->ai_screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->ai_screen_params, src->ai_screen_params, size);
	}
	if (src->polarizable_pts) {
		size = src->n_polarizable_pts * sizeof(struct polarizable_pt);
		dest->polarizable_pts = (struct polarizable_pt *)malloc(size);
		if (!dest->polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->polarizable_pts, src->polarizable_pts, size);
	}
	if (src->dynamic_polarizable_pts) {
		size = src->n_dynamic_polarizable_pts *
				sizeof(struct dynamic_polarizable_pt);
		dest->dynamic_polarizable_pts =
				(struct dynamic_polarizable_pt *)malloc(size);
		if (!dest->dynamic_polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->dynamic_polarizable_pts,
				src->dynamic_polarizable_pts, size);
	}
	if (src->lmo_centroids) {
		size = src->n_lmo * sizeof(vec_t);
		dest->lmo_centroids = (vec_t *)malloc(size);
		if (!dest->lmo_centroids)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->lmo_centroids, src->lmo_centroids, size);
	}
	if (src->xr_atoms) {
		size = src->n_xr_atoms * sizeof(struct xr_atom);
		dest->xr_atoms = (struct xr_atom *)malloc(size);
		if (!dest->xr_atoms)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_atoms, src->xr_atoms, size);

		for (size_t j = 0; j < src->n_xr_atoms; j++) {
			const struct xr_atom *at_src = src->xr_atoms + j;
			struct xr_atom *at_dest = dest->xr_atoms + j;

			size = at_src->n_shells * sizeof(struct shell);
			at_dest->shells = (struct shell *)malloc(size);
			if (!at_dest->shells)
				return EFP_RESULT_NO_MEMORY;
			memcpy(at_dest->shells, at_src->shells, size);

			for (size_t i = 0; i < at_src->n_shells; i++) {
				size = (at_src->shells[i].type == 'L' ? 3 : 2) *
					at_src->shells[i].n_funcs * sizeof(double);

				at_dest->shells[i].coef = (double *)malloc(size);
				if (!at_dest->shells[i].coef)
					return EFP_RESULT_NO_MEMORY;
				memcpy(at_dest->shells[i].coef, at_src->shells[i].coef, size);
			}
		}
	}
	if (src->xr_fock_mat) {
		size = src->n_lmo * (src->n_lmo + 1) / 2 * sizeof(double);
		dest->xr_fock_mat = (double *)malloc(size);
		if (!dest->xr_fock_mat)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_fock_mat, src->xr_fock_mat, size);
	}
	if (src->xr_wf) {
		size = src->n_lmo * src->xr_wf_size * sizeof(double);
		dest->xr_wf = (double *)malloc(size);
		if (!dest->xr_wf)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_wf, src->xr_wf, size);
	}
	if (src->xrfit) {
		size = src->n_lmo * 4 * sizeof(double);
		dest->xrfit = (double *)malloc(size);
		if (!dest->xrfit)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xrfit, src->xrfit, size);
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_opts(const struct efp_opts *opts)
{
	if (opts->enable_pbc) {
		if ((opts->terms & EFP_TERM_AI_ELEC) ||
		    (opts->terms & EFP_TERM_AI_POL) ||
		    (opts->terms & EFP_TERM_AI_DISP) ||
		    (opts->terms & EFP_TERM_AI_XR) ||
		    (opts->terms & EFP_TERM_AI_CHTR)) {
			efp_log("periodic calculations are not supported for QM/EFP");
			return EFP_RESULT_FATAL;
		}

		if (!opts->enable_cutoff) {
			efp_log("periodic calculations require interaction cutoff");
			return EFP_RESULT_FATAL;
		}
	}

	if (opts->enable_cutoff) {
		if (opts->swf_cutoff < 1.0) {
			efp_log("interaction cutoff is too small");
			return EFP_RESULT_FATAL;
		}
	}

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_frag_params(const struct efp_opts *opts, const struct frag *frag)
{
	if ((opts->terms & EFP_TERM_ELEC) || (opts->terms & EFP_TERM_AI_ELEC)) {
		if (!frag->multipole_pts) {
			efp_log("electrostatic parameters are missing");
			return EFP_RESULT_FATAL;
		}

		if (opts->elec_damp == EFP_ELEC_DAMP_SCREEN && !frag->screen_params) {
			efp_log("screening parameters are missing");
			return EFP_RESULT_FATAL;
		}
	}
	if ((opts->terms & EFP_TERM_POL) || (opts->terms & EFP_TERM_AI_POL)) {
		if (!frag->polarizable_pts || !frag->multipole_pts) {
			efp_log("polarization parameters are missing");
			return EFP_RESULT_FATAL;
		}
	}
	if ((opts->terms & EFP_TERM_DISP) || (opts->terms & EFP_TERM_AI_DISP)) {
		if (!frag->dynamic_polarizable_pts) {
			efp_log("dispersion parameters are missing");
			return EFP_RESULT_FATAL;
		}

		if (opts->disp_damp == EFP_DISP_DAMP_OVERLAP &&
		    frag->n_lmo != frag->n_dynamic_polarizable_pts) {
			efp_log("number of polarization points does not match number of LMOs");
			return EFP_RESULT_FATAL;
		}
	}
	if ((opts->terms & EFP_TERM_XR) || (opts->terms & EFP_TERM_AI_XR)) {
		if (!frag->xr_atoms ||
		    !frag->xr_fock_mat ||
		    !frag->xr_wf ||
		    !frag->lmo_centroids) {
			efp_log("exchange repulsion parameters are missing");
			return EFP_RESULT_FATAL;
		}
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_params(struct efp *efp)
{
	enum efp_result res;

	for (size_t i = 0; i < efp->n_frag; i++)
		if ((res = check_frag_params(&efp->opts, efp->frags + i)))
			return res;

	return EFP_RESULT_SUCCESS;
}

static bool
do_elec(const struct efp_opts *opts)
{
	return (opts->terms & EFP_TERM_ELEC);
}

static bool
do_disp(const struct efp_opts *opts)
{
	return (opts->terms & EFP_TERM_DISP);
}

static bool
do_xr(const struct efp_opts *opts)
{
	bool xr = (opts->terms & EFP_TERM_XR);
	bool cp = (opts->terms & EFP_TERM_ELEC) && (opts->elec_damp == EFP_ELEC_DAMP_OVERLAP);
	bool dd = (opts->terms & EFP_TERM_DISP) && (opts->disp_damp == EFP_DISP_DAMP_OVERLAP);

	return (xr || cp || dd);
}

static void
compute_two_body_range(struct efp *efp, size_t frag_from, size_t frag_to, void *data)
{
	double e_elec = 0.0, e_disp = 0.0, e_xr = 0.0, e_cp = 0.0;

	(void)data;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:e_elec,e_disp,e_xr,e_cp)
#endif
	for (size_t i = frag_from; i < frag_to; i++) {
		size_t cnt = efp_inner_count(i, efp->n_frag);

		for (size_t j = i + 1; j < i + 1 + cnt; j++) {
			size_t fr_j = j % efp->n_frag;

			if (!efp_skip_frag_pair(efp, i, fr_j)) {
				double *s;
				six_t *ds;

				size_t n_lmo_i = efp->frags[i].n_lmo;
				size_t n_lmo_j = efp->frags[fr_j].n_lmo;

				s = (double *)calloc(n_lmo_i * n_lmo_j, sizeof(double));
				ds = (six_t *)calloc(n_lmo_i * n_lmo_j, sizeof(six_t));

				if (do_xr(&efp->opts)) {
					double exr, ecp;

					efp_frag_frag_xr(efp, i, fr_j, s, ds, &exr, &ecp);
					e_xr += exr;
					e_cp += ecp;
				}

				if (do_elec(&efp->opts))
					e_elec += efp_frag_frag_elec(efp, i, fr_j);

				if (do_disp(&efp->opts))
					e_disp += efp_frag_frag_disp(efp, i, fr_j, s, ds);

				free(s);
				free(ds);
			}
		}
	}

	efp->energy.electrostatic += e_elec;
	efp->energy.dispersion += e_disp;
	efp->energy.exchange_repulsion += e_xr;
	efp->energy.charge_penetration += e_cp;
}

EFP_EXPORT enum efp_result
efp_get_energy(struct efp *efp, struct efp_energy *energy)
{
	assert(efp);
	assert(energy);

	*energy = efp->energy;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_gradient(struct efp *efp, double *grad)
{
	assert(efp);
	assert(grad);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return (EFP_RESULT_FATAL);
	}

	memcpy(grad, efp->grad, efp->n_frag * sizeof(six_t));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_point_charges(struct efp *efp, size_t n_ptc, const double *ptc, const double *xyz)
{
	assert(efp);
	efp->n_ptc = n_ptc;

	if (n_ptc == 0) {
		free(efp->ptc);
		free(efp->ptc_xyz);
		free(efp->ptc_grad);
		efp->ptc = NULL;
		efp->ptc_xyz = NULL;
		efp->ptc_grad = NULL;

		return (EFP_RESULT_SUCCESS);
	}

	assert(ptc);
	assert(xyz);

	efp->ptc = (double *)realloc(efp->ptc, n_ptc * sizeof(double));
	efp->ptc_xyz = (vec_t *)realloc(efp->ptc_xyz, n_ptc * sizeof(vec_t));
	efp->ptc_grad = (vec_t *)realloc(efp->ptc_grad, n_ptc * sizeof(vec_t));

	memcpy(efp->ptc, ptc, n_ptc * sizeof(double));
	memcpy(efp->ptc_xyz, xyz, n_ptc * sizeof(vec_t));
	memset(efp->ptc_grad, 0, n_ptc * sizeof(vec_t));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_point_charge_count(struct efp *efp, size_t *n_ptc)
{
	assert(efp);
	assert(n_ptc);

	*n_ptc = efp->n_ptc;

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_point_charge_gradient(struct efp *efp, double *grad)
{
	assert(efp);
	assert(grad);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return (EFP_RESULT_FATAL);
	}

	memcpy(grad, efp->ptc_grad, efp->n_ptc * sizeof(vec_t));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_point_charge_coordinates(struct efp *efp, double *xyz)
{
	assert(efp);
	assert(xyz);

	memcpy(xyz, efp->ptc_xyz, efp->n_ptc * sizeof(vec_t));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_point_charge_coordinates(struct efp *efp, const double *xyz)
{
	assert(efp);
	assert(xyz);

	memcpy(efp->ptc_xyz, xyz, efp->n_ptc * sizeof(vec_t));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_point_charge_values(struct efp *efp, double *ptc)
{
	assert(efp);
	assert(ptc);

	memcpy(ptc, efp->ptc, efp->n_ptc * sizeof(double));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_point_charge_values(struct efp *efp, const double *ptc)
{
	assert(efp);
	assert(ptc);

	memcpy(efp->ptc, ptc, efp->n_ptc * sizeof(double));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_coordinates(struct efp *efp, enum efp_coord_type coord_type,
			const double *coord)
{
	assert(efp);
	assert(coord);

	size_t stride;
	enum efp_result res;

	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC:
			stride = 6;
			break;
		case EFP_COORD_TYPE_POINTS:
			stride = 9;
			break;
		case EFP_COORD_TYPE_ROTMAT:
			stride = 12;
			break;
	}

	for (size_t i = 0; i < efp->n_frag; i++, coord += stride)
		if ((res = efp_set_frag_coordinates(efp, i, coord_type, coord)))
			return (res);

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_frag_coordinates(struct efp *efp, size_t frag_idx,
		enum efp_coord_type coord_type, const double *coord)
{
	struct frag *frag;

	assert(efp);
	assert(coord);
	assert(frag_idx < efp->n_frag);

	frag = efp->frags + frag_idx;

	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC:
			return (set_coord_xyzabc(frag, coord));
		case EFP_COORD_TYPE_POINTS:
			return (set_coord_points(frag, coord));
		case EFP_COORD_TYPE_ROTMAT:
			return (set_coord_rotmat(frag, coord));
	}

	/* unreachable */

	return (EFP_RESULT_FATAL);
}

EFP_EXPORT enum efp_result
efp_get_coordinates(struct efp *efp, double *xyzabc)
{
	assert(efp);
	assert(xyzabc);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		double a, b, c;
		matrix_to_euler(&frag->rotmat, &a, &b, &c);

		*xyzabc++ = frag->x;
		*xyzabc++ = frag->y;
		*xyzabc++ = frag->z;
		*xyzabc++ = a;
		*xyzabc++ = b;
		*xyzabc++ = c;
	}

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_frag_xyzabc(struct efp *efp, size_t frag_idx, double *xyzabc)
{
	struct frag *frag;
	double a, b, c;

	assert(efp);
	assert(frag_idx < efp->n_frag);
	assert(xyzabc);

	frag = efp->frags + frag_idx;
	matrix_to_euler(&frag->rotmat, &a, &b, &c);

	xyzabc[0] = frag->x;
	xyzabc[1] = frag->y;
	xyzabc[2] = frag->z;
	xyzabc[3] = a;
	xyzabc[4] = b;
	xyzabc[5] = c;

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_periodic_box(struct efp *efp, double x, double y, double z)
{
	assert(efp);

	if (x < 2.0 * efp->opts.swf_cutoff ||
	    y < 2.0 * efp->opts.swf_cutoff ||
	    z < 2.0 * efp->opts.swf_cutoff) {
		efp_log("periodic box dimensions must be at least twice the cutoff");
		return EFP_RESULT_FATAL;
	}

	efp->box.x = x;
	efp->box.y = y;
	efp->box.z = z;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_stress_tensor(struct efp *efp, double *stress)
{
	assert(efp);
	assert(stress);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return (EFP_RESULT_FATAL);
	}

	*(mat_t *)stress = efp->stress;

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_prepare(struct efp *efp)
{
	assert(efp);

	efp->n_polarizable_pts = 0;

	for (size_t i = 0; i < efp->n_frag; i++) {
		efp->frags[i].polarizable_offset = efp->n_polarizable_pts;
		efp->n_polarizable_pts += efp->frags[i].n_polarizable_pts;
	}

	efp->indip = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));
	efp->indipconj = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));
	efp->grad = (six_t *)calloc(efp->n_frag, sizeof(six_t));

	if ((efp->skiplist = efp_bvec_create(efp->n_frag * efp->n_frag)) == NULL)
		return (EFP_RESULT_NO_MEMORY);

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_orbital_energies(struct efp *efp, size_t n_core, size_t n_act,
				size_t n_vir, const double *oe)
{
	size_t size;

	assert(efp);
	assert(oe);

	efp->n_ai_core = n_core;
	efp->n_ai_act = n_act;
	efp->n_ai_vir = n_vir;

	size = (n_core + n_act + n_vir) * sizeof(double);

	efp->ai_orbital_energies = (double *)realloc(efp->ai_orbital_energies, size);
	memcpy(efp->ai_orbital_energies, oe, size);

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_set_dipole_integrals(struct efp *efp, size_t n_core, size_t n_act,
				size_t n_vir, const double *dipint)
{
	size_t size;

	assert(efp);
	assert(dipint);

	efp->n_ai_core = n_core;
	efp->n_ai_act = n_act;
	efp->n_ai_vir = n_vir;

	size = 3 * (n_core + n_act + n_vir) * (n_core + n_act + n_vir) * sizeof(double);

	efp->ai_dipole_integrals = (double *)realloc(efp->ai_dipole_integrals, size);
	memcpy(efp->ai_dipole_integrals, dipint, size);

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_wavefunction_dependent_energy(struct efp *efp, double *energy)
{
	assert(efp);
	assert(energy);

	if (!(efp->opts.terms & EFP_TERM_POL) && !(efp->opts.terms & EFP_TERM_AI_POL)) {
		*energy = 0.0;
		return EFP_RESULT_SUCCESS;
	}

	return efp_compute_pol_energy(efp, energy);
}

EFP_EXPORT enum efp_result
efp_compute(struct efp *efp, int do_gradient)
{
	enum efp_result res;

	assert(efp);

	if (efp->grad == NULL) {
		efp_log("call efp_prepare after all fragments are added");
		return (EFP_RESULT_FATAL);
	}

	efp->do_gradient = do_gradient;

	if ((res = check_params(efp)))
		return (res);

	memset(&efp->energy, 0, sizeof(struct efp_energy));
	memset(&efp->stress, 0, sizeof(mat_t));
	memset(efp->grad, 0, efp->n_frag * sizeof(six_t));
	memset(efp->ptc_grad, 0, efp->n_ptc * sizeof(vec_t));

	efp_balance_work(efp, compute_two_body_range, NULL);

	if ((res = efp_compute_pol(efp)))
		return res;

	if ((res = efp_compute_ai_elec(efp)))
		return res;

	if ((res = efp_compute_ai_disp(efp)))
		return res;

#ifdef WITH_MPI
	efp_allreduce(&efp->energy.electrostatic, 1);
	efp_allreduce(&efp->energy.dispersion, 1);
	efp_allreduce(&efp->energy.exchange_repulsion, 1);
	efp_allreduce(&efp->energy.charge_penetration, 1);

	if (efp->do_gradient) {
		efp_allreduce((double *)efp->grad, 6 * efp->n_frag);
		efp_allreduce((double *)efp->ptc_grad, 3 * efp->n_ptc);
		efp_allreduce((double *)&efp->stress, 9);
	}
#endif
	efp->energy.total = efp->energy.electrostatic +
			    efp->energy.charge_penetration +
			    efp->energy.electrostatic_point_charges +
			    efp->energy.polarization +
			    efp->energy.dispersion +
			    efp->energy.ai_dispersion +
			    efp->energy.exchange_repulsion;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_charge(struct efp *efp, size_t frag_idx, double *charge)
{
	assert(efp);
	assert(charge);
	assert(frag_idx < efp->n_frag);

	struct frag *frag = efp->frags + frag_idx;
	*charge = 0.0;

	for (size_t i = 0; i < frag->n_atoms; i++)
		*charge += frag->atoms[i].znuc;

	for (size_t i = 0; i < frag->n_multipole_pts; i++)
		*charge += frag->multipole_pts[i].monopole;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_multiplicity(struct efp *efp, size_t frag_idx, int *mult)
{
	assert(efp);
	assert(mult);
	assert(frag_idx < efp->n_frag);

	*mult = efp->frags[frag_idx].multiplicity;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_count(struct efp *efp, size_t *n_mult)
{
	size_t sum = 0;

	assert(efp);
	assert(n_mult);

	for (size_t i = 0; i < efp->n_frag; i++)
		sum += efp->frags[i].n_multipole_pts;

	*n_mult = sum;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_coordinates(struct efp *efp, double *xyz)
{
	assert(efp);
	assert(xyz);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_multipole_pts; j++) {
			*xyz++ = frag->multipole_pts[j].x;
			*xyz++ = frag->multipole_pts[j].y;
			*xyz++ = frag->multipole_pts[j].z;
		}
	}

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_values(struct efp *efp, double *mult)
{
	assert(efp);
	assert(mult);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_multipole_pts; j++) {
			struct multipole_pt *pt = frag->multipole_pts + j;

			*mult++ = pt->monopole;

			*mult++ = pt->dipole.x;
			*mult++ = pt->dipole.y;
			*mult++ = pt->dipole.z;

			for (size_t t = 0; t < 6; t++)
				*mult++ = pt->quadrupole[t];

			for (size_t t = 0; t < 10; t++)
				*mult++ = pt->octupole[t];
		}
	}

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_count(struct efp *efp, size_t *n_dip)
{
	size_t sum = 0;

	assert(efp);
	assert(n_dip);

	for (size_t i = 0; i < efp->n_frag; i++)
		sum += efp->frags[i].n_polarizable_pts;

	*n_dip = sum;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_coordinates(struct efp *efp, double *xyz)
{
	assert(efp);
	assert(xyz);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			*xyz++ = pt->x;
			*xyz++ = pt->y;
			*xyz++ = pt->z;
		}
	}

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_values(struct efp *efp, double *dip)
{
	assert(efp);
	assert(dip);

	memcpy(dip, efp->indip, efp->n_polarizable_pts * sizeof(vec_t));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_conj_values(struct efp *efp, double *dip)
{
	assert(efp);
	assert(dip);

	memcpy(dip, efp->indipconj, efp->n_polarizable_pts * sizeof(vec_t));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_lmo_count(struct efp *efp, size_t frag_idx, size_t *n_lmo)
{
	assert(efp != NULL);
	assert(frag_idx < efp->n_frag);
	assert(n_lmo != NULL);

	*n_lmo = efp->frags[frag_idx].n_lmo;

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_lmo_coordinates(struct efp *efp, size_t frag_idx, double *xyz)
{
	struct frag *frag;

	assert(efp != NULL);
	assert(frag_idx < efp->n_frag);
	assert(xyz != NULL);

	frag = efp->frags + frag_idx;

	if (frag->lmo_centroids == NULL) {
		efp_log("no LMO centroids for fragment %s", frag->name);
		return (EFP_RESULT_FATAL);
	}

	memcpy(xyz, frag->lmo_centroids, frag->n_lmo * sizeof(vec_t));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT enum efp_result
efp_get_xrfit(struct efp *efp, size_t frag_idx, double *xrfit)
{
	struct frag *frag;

	assert(efp != NULL);
	assert(frag_idx < efp->n_frag);
	assert(xrfit != NULL);

	frag = efp->frags + frag_idx;

	if (frag->xrfit == NULL) {
		efp_log("no XRFIT parameters for fragment %s", frag->name);
		return (EFP_RESULT_FATAL);
	}

	memcpy(xrfit, frag->xrfit, frag->n_lmo * 4 * sizeof(double));

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT void
efp_shutdown(struct efp *efp)
{
	if (!efp)
		return;

	for (size_t i = 0; i < efp->n_frag; i++)
		free_frag(efp->frags + i);

	for (size_t i = 0; i < efp->n_lib; i++) {
		free_frag(efp->lib[i]);
		free(efp->lib[i]);
	}

	free(efp->frags);
	free(efp->lib);
	free(efp->grad);
	free(efp->ptc);
	free(efp->ptc_xyz);
	free(efp->ptc_grad);
	free(efp->indip);
	free(efp->indipconj);
	free(efp->ai_orbital_energies);
	free(efp->ai_dipole_integrals);
	efp_bvec_free(efp->skiplist);
	free(efp);
}

EFP_EXPORT enum efp_result
efp_set_opts(struct efp *efp, const struct efp_opts *opts)
{
	enum efp_result res;

	assert(efp);
	assert(opts);

	if ((res = check_opts(opts)))
		return res;

	efp->opts = *opts;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_opts(struct efp *efp, struct efp_opts *opts)
{
	assert(efp);
	assert(opts);

	*opts = efp->opts;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT void
efp_opts_default(struct efp_opts *opts)
{
	assert(opts);

	memset(opts, 0, sizeof(struct efp_opts));

	opts->terms = EFP_TERM_ELEC | EFP_TERM_POL | EFP_TERM_DISP |
		EFP_TERM_XR | EFP_TERM_AI_ELEC | EFP_TERM_AI_POL;
}

EFP_EXPORT void
efp_set_error_log(void (*cb)(const char *))
{
	efp_set_log_cb(cb);
}

EFP_EXPORT enum efp_result
efp_add_fragment(struct efp *efp, const char *name)
{
	assert(efp);
	assert(name);

	enum efp_result res;
	const struct frag *lib = efp_find_lib(efp, name);

	if (!lib)
		return EFP_RESULT_UNKNOWN_FRAGMENT;

	efp->n_frag++;
	efp->frags = (struct frag *)realloc(efp->frags,
		efp->n_frag * sizeof(struct frag));

	if (!efp->frags)
		return EFP_RESULT_NO_MEMORY;

	struct frag *frag = efp->frags + efp->n_frag - 1;

	if ((res = copy_frag(frag, lib)))
		return res;

	for (size_t a = 0; a < 3; a++) {
		size_t size = frag->xr_wf_size * frag->n_lmo;

		frag->xr_wf_deriv[a] = (double *)calloc(size, sizeof(double));

		if (!frag->xr_wf_deriv[a])
			return EFP_RESULT_NO_MEMORY;
	}

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_skip_fragments(struct efp *efp, size_t i, size_t j, int value)
{
	assert(efp);
	assert(efp->skiplist); /* call efp_prepare first */
	assert(i < efp->n_frag);
	assert(j < efp->n_frag);

	efp_bvec_set_value(efp->skiplist, i * efp->n_frag + j, value);
	efp_bvec_set_value(efp->skiplist, j * efp->n_frag + i, value);

	return (EFP_RESULT_SUCCESS);
}

EFP_EXPORT struct efp *
efp_create(void)
{
	struct efp *efp = (struct efp *)calloc(1, sizeof(struct efp));

	if (!efp)
		return (NULL);

	efp_opts_default(&efp->opts);

	return (efp);
}

EFP_EXPORT enum efp_result
efp_set_electron_density_field_fn(struct efp *efp, efp_electron_density_field_fn fn)
{
	assert(efp);

	efp->get_electron_density_field = fn;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_electron_density_field_user_data(struct efp *efp, void *user_data)
{
	assert(efp);

	efp->get_electron_density_field_user_data = user_data;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_count(struct efp *efp, size_t *n_frag)
{
	assert(efp);
	assert(n_frag);

	*n_frag = efp->n_frag;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_name(struct efp *efp, size_t frag_idx, size_t size, char *frag_name)
{
	assert(efp);
	assert(frag_name);
	assert(frag_idx < efp->n_frag);

	strncpy(frag_name, efp->frags[frag_idx].name, size);

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_mass(struct efp *efp, size_t frag_idx, double *mass_out)
{
	assert(efp);
	assert(mass_out);
	assert(frag_idx < efp->n_frag);

	const struct frag *frag = efp->frags + frag_idx;
	double mass = 0.0;

	for (size_t i = 0; i < frag->n_atoms; i++)
		mass += frag->atoms[i].mass;

	*mass_out = mass;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_inertia(struct efp *efp, size_t frag_idx, double *inertia_out)
{
	assert(efp);
	assert(inertia_out);
	assert(frag_idx < efp->n_frag);

	/* center of mass is in origin and axes are principal axes of inertia */

	const struct frag *frag = efp->frags[frag_idx].lib;
	vec_t inertia = vec_zero;

	for (size_t i = 0; i < frag->n_atoms; i++) {
		const struct efp_atom *atom = frag->atoms + i;

		inertia.x += atom->mass * (atom->y * atom->y + atom->z * atom->z);
		inertia.y += atom->mass * (atom->x * atom->x + atom->z * atom->z);
		inertia.z += atom->mass * (atom->x * atom->x + atom->y * atom->y);
	}

	inertia_out[0] = inertia.x;
	inertia_out[1] = inertia.y;
	inertia_out[2] = inertia.z;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_atom_count(struct efp *efp, size_t frag_idx, size_t *n_atoms)
{
	assert(efp);
	assert(n_atoms);
	assert(frag_idx < efp->n_frag);

	*n_atoms = efp->frags[frag_idx].n_atoms;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_atoms(struct efp *efp, size_t frag_idx, size_t size, struct efp_atom *atoms)
{
	assert(efp);
	assert(atoms);
	assert(frag_idx < efp->n_frag);
	assert(size >= efp->frags[frag_idx].n_atoms);

	struct frag *frag = efp->frags + frag_idx;

	memcpy(atoms, frag->atoms, frag->n_atoms * sizeof(struct efp_atom));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT void
efp_torque_to_derivative(const double *euler, const double *torque, double *deriv)
{
	assert(euler);
	assert(torque);
	assert(deriv);

	double tx = torque[0];
	double ty = torque[1];
	double tz = torque[2];

	double sina = sin(euler[0]);
	double cosa = cos(euler[0]);
	double sinb = sin(euler[1]);
	double cosb = cos(euler[1]);

	deriv[0] = tz;
	deriv[1] = cosa * tx + sina * ty;
	deriv[2] = sinb * sina * tx - sinb * cosa * ty + cosb * tz;
}

EFP_EXPORT const char *
efp_banner(void)
{
	static const char banner[] =
		"LIBEFP ver. " LIBEFP_VERSION_STRING "\n"
		"Copyright (c) 2012-2014 Ilya Kaliman\n"
		"\n"
		"Journal Reference:\n"
		"    Kaliman and Slipchenko, JCC 2013.\n"
		"    DOI: http://dx.doi.org/10.1002/jcc.23375\n"
		"\n"
		"Project web site: http://www.libefp.org/\n";

	return banner;
}

EFP_EXPORT const char *
efp_result_to_string(enum efp_result res)
{
	switch (res) {
	case EFP_RESULT_SUCCESS:
		return "Operation was successful.";
	case EFP_RESULT_FATAL:
		return "Fatal error has occurred.";
	case EFP_RESULT_NO_MEMORY:
		return "Insufficient memory.";
	case EFP_RESULT_FILE_NOT_FOUND:
		return "File not found.";
	case EFP_RESULT_SYNTAX_ERROR:
		return "Syntax error.";
	case EFP_RESULT_UNKNOWN_FRAGMENT:
		return "Unknown EFP fragment.";
	case EFP_RESULT_POL_NOT_CONVERGED:
		return "Polarization SCF procedure did not converge.";
	}
	assert(0);
}
