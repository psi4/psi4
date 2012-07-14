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

#include "efp_private.h"
#include "elec.h"
#include "parse.h"

static inline int
initialized(struct efp *efp)
{
	return efp && efp->magic == EFP_INIT_MAGIC;
}

EFP_EXPORT enum efp_result
efp_get_energy(struct efp *efp, struct efp_energy *energy)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!energy)
		return EFP_RESULT_ARGUMENT_NULL;

	memcpy(energy, &efp->energy, sizeof(struct efp_energy));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_gradient(struct efp *efp, int size, double *grad)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!grad)
		return EFP_RESULT_ARGUMENT_NULL;

	if (!efp->do_gradient)
		return EFP_RESULT_GRADIENT_NOT_REQUESTED;

	if (size < 6 * efp->n_frag)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	for (int i = 0; i < efp->n_frag; i++) {
		*grad++ = efp->frags[i].force.x;
		*grad++ = efp->frags[i].force.y;
		*grad++ = efp->frags[i].force.z;
		*grad++ = efp->frags[i].torque.x;
		*grad++ = efp->frags[i].torque.y;
		*grad++ = efp->frags[i].torque.z;
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_qm_gradient(struct efp *efp, int size, double *grad)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!grad)
		return EFP_RESULT_ARGUMENT_NULL;

	if (!efp->do_gradient)
		return EFP_RESULT_GRADIENT_NOT_REQUESTED;

	if (size < 3 * efp->qm.n_atoms)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	memcpy(grad, efp->qm.grad, 3 * efp->qm.n_atoms * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_qm_atoms(struct efp *efp, int n_atoms,
		 const double *znuc, const double *xyz)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!znuc || !xyz)
		return EFP_RESULT_ARGUMENT_NULL;

	if (n_atoms != efp->qm.n_atoms) {
		efp->qm.znuc =
			realloc(efp->qm.znuc, n_atoms * sizeof(double));
		efp->qm.xyz =
			realloc(efp->qm.xyz, 3 * n_atoms * sizeof(double));
		efp->qm.grad =
			realloc(efp->qm.grad, 3 * n_atoms * sizeof(double));
	}

	memcpy(efp->qm.znuc, znuc, n_atoms * sizeof(double));
	memcpy(efp->qm.xyz, xyz, 3 * n_atoms * sizeof(double));

	return EFP_RESULT_SUCCESS;
}

static void
euler_to_matrix(double a, double b, double c, mat_t *out)
{
	double sina = sin(a), cosa = cos(a);
	double sinb = sin(b), cosb = cos(b);
	double sinc = sin(c), cosc = cos(c);
	out->xx =  cosa * cosc - sina * cosb * sinc;
	out->xy = -cosa * sinc - sina * cosb * cosc;
	out->xz =  sinb * sina;
	out->yx =  sina * cosc + cosa * cosb * sinc;
	out->yy = -sina * sinc + cosa * cosb * cosc;
	out->yz = -sinb * cosa;
	out->zx =  sinb * sinc;
	out->zy =  sinb * cosc;
	out->zz =  cosb;
}

static void
matrix_to_euler(const mat_t *rotmat, double *ea, double *eb, double *ec)
{
	double a, b, c, sinb;

	b = acos(rotmat->zz);
	sinb = sqrt(1.0 - rotmat->zz * rotmat->zz);

	if (fabs(sinb) < 1.0e-6) {
		a = atan2(-rotmat->xy, rotmat->xx);
		c = 0.0;
	}
	else {
		a = atan2(rotmat->xz, -rotmat->yz);
		c = atan2(rotmat->zx, rotmat->zy);
	}

	*ea = a, *eb = b, *ec = c;
}

static void
update_fragment(struct frag *frag)
{
	/* update atoms */
	for (int i = 0; i < frag->n_atoms; i++)
		move_pt(VEC(frag->x), &frag->rotmat, VEC(frag->lib->x),
			VEC(frag->lib->atoms[i].x), VEC(frag->atoms[i].x));

	efp_update_elec(frag);
	efp_update_pol(frag);
	efp_update_disp(frag);
	efp_update_xr(frag);
}

EFP_EXPORT enum efp_result
efp_set_coordinates(struct efp *efp, const double *xyzabc)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!xyzabc)
		return EFP_RESULT_ARGUMENT_NULL;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		frag->x = *xyzabc++;
		frag->y = *xyzabc++;
		frag->z = *xyzabc++;

		double a = *xyzabc++;
		double b = *xyzabc++;
		double c = *xyzabc++;

		euler_to_matrix(a, b, c, &frag->rotmat);

		update_fragment(frag);
	}
	return EFP_RESULT_SUCCESS;
}

static void
points_to_matrix(const double *pts, mat_t *rotmat)
{
	double (*rm)[3] = (double (*)[3])rotmat;

	const double *p1 = pts + 0;
	const double *p2 = pts + 3;
	const double *p3 = pts + 6;

	double t1norm = 0.0;
	double t2norm = 0.0;

	for (int i = 0; i < 3; i++) {
		rm[i][0] = p2[i] - p1[i];
		t1norm += rm[i][0] * rm[i][0];
		rm[i][1] = p3[i] - p1[i];
		t2norm += rm[i][1] * rm[i][1];
	}

	t1norm = 1.0 / sqrt(t1norm);
	t2norm = 1.0 / sqrt(t2norm);

	for (int i = 0; i < 3; i++) {
		rm[i][0] *= t1norm;
		rm[i][1] *= t2norm;
	}

	double dot = rm[0][0] * rm[0][1] +
		     rm[1][0] * rm[1][1] +
		     rm[2][0] * rm[2][1];

	rm[0][1] -= dot * rm[0][0];
	rm[1][1] -= dot * rm[1][0];
	rm[2][1] -= dot * rm[2][0];

	rm[0][2] = rm[1][0] * rm[2][1] - rm[2][0] * rm[1][1];
	rm[1][2] = rm[2][0] * rm[0][1] - rm[0][0] * rm[2][1];
	rm[2][2] = rm[0][0] * rm[1][1] - rm[1][0] * rm[0][1];

	for (int j = 0; j < 3; j++) {
		double vecsq = 0.0;
		for (int i = 0; i < 3; i++)
			vecsq += rm[i][j] * rm[i][j];
		vecsq = sqrt(vecsq);
		for (int i = 0; i < 3; i++)
			rm[i][j] /= vecsq;
	}
}

EFP_EXPORT enum efp_result
efp_set_coordinates_2(struct efp *efp, const double *pts)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!pts)
		return EFP_RESULT_ARGUMENT_NULL;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		const double *pt = pts + 9 * i;

		points_to_matrix(pt, &frag->rotmat);

		vec_t p1;
		mat_vec(&frag->rotmat, VEC(frag->lib->atoms[0].x), &p1);

		/* center of mass */
		frag->x = pt[0] - p1.x;
		frag->y = pt[1] - p1.y;
		frag->z = pt[2] - p1.z;

		update_fragment(frag);
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_coordinates(struct efp *efp, int size, double *xyzabc)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!xyzabc)
		return EFP_RESULT_ARGUMENT_NULL;

	if (size < 6 * efp->n_frag)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	for (int i = 0; i < efp->n_frag; i++) {
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
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_scf_update(struct efp *efp, double *energy)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!energy)
		return EFP_RESULT_ARGUMENT_NULL;

	*energy = efp_compute_pol_energy(efp);
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_compute(struct efp *efp, int do_gradient)
{
	typedef enum efp_result (*term_fn)(struct efp *);

	static const term_fn term_list[] = {
		efp_compute_elec,
		efp_compute_pol,
		efp_compute_disp,
		efp_compute_xr,
		efp_compute_chtr,
		efp_compute_ai_elec,
		efp_compute_ai_disp,
		efp_compute_ai_xr,
		efp_compute_ai_chtr
	};

	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	efp->do_gradient = do_gradient;

	if (do_gradient) {
		for (int i = 0; i < efp->n_frag; i++) {
			vec_zero(&efp->frags[i].force);
			vec_zero(&efp->frags[i].torque);
		}
		memset(efp->qm.grad, 0, 3 * efp->qm.n_atoms * sizeof(double));
	}

	enum efp_result res;

	for (size_t i = 0; i < ARRAY_SIZE(term_list); i++)
		if ((res = term_list[i](efp)))
			return res;

	efp->energy.total = efp->energy.electrostatic +
			    efp->energy.charge_penetration +
			    efp->energy.polarization +
			    efp->energy.dispersion +
			    efp->energy.exchange_repulsion +
			    efp->energy.charge_transfer +
			    efp->energy.ai_electrostatic +
			    efp->energy.ai_dispersion +
			    efp->energy.ai_exchange_repulsion +
			    efp->energy.ai_charge_transfer;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_count(struct efp *efp, int *n_mult)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!n_mult)
		return EFP_RESULT_ARGUMENT_NULL;

	int n_charge = 0;
	int n_dipole = 0;
	int n_quadrupole = 0;
	int n_octupole = 0;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		n_charge += frag->n_atoms;
		n_charge += frag->n_multipole_pts;

		n_dipole += frag->n_polarizable_pts;
		n_dipole += frag->n_multipole_pts;

		n_quadrupole += frag->n_multipole_pts;

		n_octupole += frag->n_multipole_pts;
	}

	n_mult[0] = n_charge;
	n_mult[1] = n_dipole;
	n_mult[2] = n_quadrupole;
	n_mult[3] = n_octupole;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipoles(struct efp *efp, double **xyz, double **z)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!xyz || !z)
		return EFP_RESULT_ARGUMENT_NULL;

	double *xyz_c = xyz[0];
	double *xyz_d = xyz[1];
	double *xyz_q = xyz[2];
	double *xyz_o = xyz[3];

	double *z_c = z[0];
	double *z_d = z[1];
	double *z_q = z[2];
	double *z_o = z[3];

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		/* atom charges */
		for (int j = 0; j < frag->n_atoms; j++) {
			struct efp_atom *at = frag->atoms + j;

			*xyz_c++ = at->x;
			*xyz_c++ = at->y;
			*xyz_c++ = at->z;

			*z_c++ = at->znuc;
		}

		/* induced dipoles */
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			*xyz_d++ = pt->x;
			*xyz_d++ = pt->y;
			*xyz_d++ = pt->z;

			*z_d++ = 0.5 * (pt->induced_dipole.x +
					pt->induced_dipole_conj.x);
			*z_d++ = 0.5 * (pt->induced_dipole.y +
					pt->induced_dipole_conj.y);
			*z_d++ = 0.5 * (pt->induced_dipole.z +
					pt->induced_dipole_conj.z);
		}
	}

	/* multipoles from electrostatics */
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (int j = 0; j < frag->n_multipole_pts; j++) {
			struct multipole_pt *pt = frag->multipole_pts + j;

			*xyz_c++ = *xyz_d++ = *xyz_q++ = *xyz_o++ = pt->x;
			*xyz_c++ = *xyz_d++ = *xyz_q++ = *xyz_o++ = pt->y;
			*xyz_c++ = *xyz_d++ = *xyz_q++ = *xyz_o++ = pt->z;

			*z_c++ = pt->monopole;

			*z_d++ = pt->dipole.x;
			*z_d++ = pt->dipole.y;
			*z_d++ = pt->dipole.z;

			for (int t = 0; t < 6; t++)
				*z_q++ = pt->quadrupole[t];

			for (int t = 0; t < 10; t++)
				*z_o++ = pt->octupole[t];
		}
	}

	return EFP_RESULT_SUCCESS;
}

static void
free_frag(struct frag *frag)
{
	if (!frag)
		return;

	free(frag->name);
	free(frag->atoms);
	free(frag->multipole_pts);
	free(frag->polarizable_pts);
	free(frag->dynamic_polarizable_pts);
	free(frag->lmo_centroids);
	free(frag->xr_fock_mat);
	free(frag->xr_wf);
	free(frag->screen_params);
	free(frag->ai_screen_params);

	for (int i = 0; i < frag->n_xr_shells; i++)
		free(frag->xr_shells[i].coef);

	free(frag->xr_shells);

	/* don't do free(frag) here */
}

EFP_EXPORT void
efp_shutdown(struct efp *efp)
{
	if (!efp)
		return;

	for (int i = 0; i < efp->n_frag; i++)
		free_frag(efp->frags + i);

	for (int i = 0; i < efp->n_lib; i++)
		free_frag(efp->lib + i);

	free(efp->frags);
	free(efp->lib);
	free(efp->disp_damp_overlap_offset);
	free(efp->disp_damp_overlap);
	free(efp->qm.znuc);
	free(efp->qm.xyz);
	free(efp->qm.grad);
	free(efp);
}

static enum efp_result
copy_frag(struct frag *dest, const struct frag *src)
{
	size_t size;
	memcpy(dest, src, sizeof(struct frag));

	if (src->name) {
		dest->name = strdup(src->name);
		if (!dest->name)
			return EFP_RESULT_NO_MEMORY;
	}
	if (src->atoms) {
		size = src->n_atoms * sizeof(struct efp_atom);
		dest->atoms = malloc(size);
		if (!dest->atoms)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->atoms, src->atoms, size);
	}
	if (src->multipole_pts) {
		size = src->n_multipole_pts * sizeof(struct multipole_pt);
		dest->multipole_pts = malloc(size);
		if (!dest->multipole_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->multipole_pts, src->multipole_pts, size);
	}
	if (src->screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->screen_params = malloc(size);
		if (!dest->screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->screen_params, src->screen_params, size);
	}
	if (src->ai_screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->ai_screen_params = malloc(size);
		if (!dest->ai_screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->ai_screen_params, src->ai_screen_params, size);
	}
	if (src->polarizable_pts) {
		size = src->n_polarizable_pts * sizeof(struct polarizable_pt);
		dest->polarizable_pts = malloc(size);
		if (!dest->polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->polarizable_pts, src->polarizable_pts, size);
	}
	if (src->dynamic_polarizable_pts) {
		size = src->n_dynamic_polarizable_pts *
				sizeof(struct dynamic_polarizable_pt);
		dest->dynamic_polarizable_pts = malloc(size);
		if (!dest->dynamic_polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->dynamic_polarizable_pts,
				src->dynamic_polarizable_pts, size);
	}
	if (src->lmo_centroids) {
		size = src->n_lmo * sizeof(vec_t);
		dest->lmo_centroids = malloc(size);
		if (!dest->lmo_centroids)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->lmo_centroids, src->lmo_centroids, size);
	}
	if (src->xr_shells) {
		size = src->n_xr_shells * sizeof(struct shell);
		dest->xr_shells = malloc(size);
		if (!dest->xr_shells)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_shells, src->xr_shells, size);

		for (int i = 0; i < src->n_xr_shells; i++) {
			size = (src->xr_shells[i].type == 'L' ? 3 : 2) *
				src->xr_shells[i].n_funcs * sizeof(double);

			dest->xr_shells[i].coef = malloc(size);
			if (!dest->xr_shells[i].coef)
				return EFP_RESULT_NO_MEMORY;
			memcpy(dest->xr_shells[i].coef,
					src->xr_shells[i].coef, size);
		}
	}
	if (src->xr_fock_mat) {
		size = src->n_lmo * (src->n_lmo + 1) / 2 * sizeof(double);
		dest->xr_fock_mat = malloc(size);
		if (!dest->xr_fock_mat)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_fock_mat, src->xr_fock_mat, size);
	}
	if (src->xr_wf) {
		size = src->n_lmo * src->xr_wf_size * sizeof(double);
		dest->xr_wf = malloc(size);
		if (!dest->xr_wf)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_wf, src->xr_wf, size);
	}
	return EFP_RESULT_SUCCESS;
}

static struct frag *
find_frag_in_library(struct efp *efp, const char *name)
{
	for (int i = 0; i < efp->n_lib; i++)
		if (!strcasecmp(efp->lib[i].name, name))
			return efp->lib + i;

	return NULL;
}

EFP_EXPORT void
efp_opts_default(struct efp_opts *opts)
{
	if (!opts)
		return;

	memset(opts, 0, sizeof(struct efp_opts));

	opts->terms = EFP_TERM_ELEC | EFP_TERM_POL | EFP_TERM_DISP |
		EFP_TERM_XR | EFP_TERM_AI_ELEC | EFP_TERM_AI_POL;
}

static enum efp_result
check_opts(const struct efp_opts *opts)
{
	unsigned terms = opts->terms;

	if (((terms & EFP_TERM_AI_ELEC) && !(terms & EFP_TERM_ELEC)) ||
	    ((terms & EFP_TERM_AI_POL) && !(terms & EFP_TERM_POL)) ||
	    ((terms & EFP_TERM_POL) && !(terms & EFP_TERM_ELEC)) ||
	    ((terms & EFP_TERM_AI_DISP) && !(terms & EFP_TERM_DISP)) ||
	    ((terms & EFP_TERM_AI_XR) && !(terms & EFP_TERM_XR)) ||
	    ((terms & EFP_TERM_AI_CHTR) && !(terms & EFP_TERM_CHTR)))
		return EFP_RESULT_INCONSISTENT_TERMS;

	if (terms & EFP_TERM_ELEC) {
		if (opts->elec_damp == EFP_ELEC_DAMP_OVERLAP &&
				!(terms & EFP_TERM_XR))
			return EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED;
	}
	if (terms & EFP_TERM_DISP) {
		if (opts->disp_damp == EFP_DISP_DAMP_OVERLAP &&
				!(terms & EFP_TERM_XR))
			return EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED;
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_params(struct efp *efp)
{
	if (efp->opts.terms & EFP_TERM_ELEC) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].multipole_pts)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	if (efp->opts.terms & EFP_TERM_POL) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].polarizable_pts)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	if (efp->opts.terms & EFP_TERM_DISP) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].dynamic_polarizable_pts)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	if (efp->opts.terms & EFP_TERM_XR) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].xr_fock_mat ||
			    !efp->frags[i].xr_wf ||
			    !efp->frags[i].lmo_centroids)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
setup_disp(struct efp *efp)
{
	if (efp->opts.disp_damp != EFP_DISP_DAMP_OVERLAP)
		return EFP_RESULT_SUCCESS;

	int n_disp = 0;
	for (int i = 0; i < efp->n_frag; i++)
		n_disp += efp->frags[i].n_dynamic_polarizable_pts;

	if (n_disp == 0)
		return EFP_RESULT_SUCCESS;

	efp->disp_damp_overlap_offset = malloc((efp->n_frag + 1) * sizeof(int));
	if (!efp->disp_damp_overlap_offset)
		return EFP_RESULT_NO_MEMORY;

	efp->disp_damp_overlap_offset[0] = 0;
	for (int i = 1; i <= efp->n_frag; i++)
		efp->disp_damp_overlap_offset[i] =
			efp->disp_damp_overlap_offset[i - 1] +
			efp->frags[i - 1].n_dynamic_polarizable_pts;

	/* XXX - this needs a lot of memory */
	efp->disp_damp_overlap = malloc(n_disp * n_disp * sizeof(double));
	if (!efp->disp_damp_overlap)
		return EFP_RESULT_NO_MEMORY;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_init(struct efp **out,
	 const struct efp_opts *opts,
	 const struct efp_callbacks *callbacks,
	 const char **potential_file_list,
	 const char **frag_name_list)
{
	if (!out || !opts || !potential_file_list || !frag_name_list)
		return EFP_RESULT_ARGUMENT_NULL;

	*out = calloc(1, sizeof(struct efp));
	if (!*out)
		return EFP_RESULT_NO_MEMORY;

	enum efp_result res;
	struct efp *efp = *out;

	if ((res = check_opts(opts)))
		return res;

	memcpy(&efp->opts, opts, sizeof(struct efp_opts));

	if (callbacks)
		memcpy(&efp->callbacks, callbacks,
				sizeof(struct efp_callbacks));

	if ((res = efp_read_potential(efp, potential_file_list)))
		return res;

	while (frag_name_list[efp->n_frag])
		efp->n_frag++;

	efp->frags = calloc(efp->n_frag, sizeof(struct frag));
	if (!efp->frags)
		return EFP_RESULT_NO_MEMORY;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = find_frag_in_library(efp,
							 frag_name_list[i]);
		if (!frag)
			return EFP_RESULT_UNKNOWN_FRAGMENT;

		if ((res = copy_frag(efp->frags + i, frag)))
			return res;
	}

	if ((res = setup_disp(efp)))
		return res;

	if ((res = check_params(efp)))
		return res;

	efp->magic = EFP_INIT_MAGIC;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_count(struct efp *efp, int *n_frag)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!n_frag)
		return EFP_RESULT_ARGUMENT_NULL;

	*n_frag = efp->n_frag;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_name(struct efp *efp, int frag_idx, int size, char *frag_name)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!frag_name)
		return EFP_RESULT_ARGUMENT_NULL;

	if (frag_idx < 0 || frag_idx >= efp->n_frag)
		return EFP_RESULT_INDEX_OUT_OF_RANGE;

	if ((unsigned)size <= strlen(efp->frags[frag_idx].name))
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	strcpy(frag_name, efp->frags[frag_idx].name);
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_atom_count(struct efp *efp, int frag_idx, int *n_atoms)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!n_atoms)
		return EFP_RESULT_ARGUMENT_NULL;

	if (frag_idx < 0 || frag_idx >= efp->n_frag)
		return EFP_RESULT_INDEX_OUT_OF_RANGE;

	*n_atoms = efp->frags[frag_idx].n_atoms;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_atoms(struct efp *efp, int frag_idx,
		   int size, struct efp_atom *atoms)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!atoms)
		return EFP_RESULT_ARGUMENT_NULL;

	if (frag_idx < 0 || frag_idx >= efp->n_frag)
		return EFP_RESULT_INDEX_OUT_OF_RANGE;

	struct frag *frag = efp->frags + frag_idx;

	if (size < frag->n_atoms)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	memcpy(atoms, frag->atoms, frag->n_atoms * sizeof(struct efp_atom));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT const char *
efp_banner(void)
{
	static const char banner[] =
		"libefp\n"
		"The Effective Fragment Potential method implementation\n"
		"Copyright (c) 2012 Ilya Kaliman\n"
		"See LICENSE file for licensing terms\n"
		"Project web page http://libefp.github.com/\n";

	return banner;
}

EFP_EXPORT const char *
efp_result_to_string(enum efp_result res)
{
	switch (res) {
	case EFP_RESULT_SUCCESS:
return "no error";
	case EFP_RESULT_NO_MEMORY:
return "out of memory";
	case EFP_RESULT_NOT_IMPLEMENTED:
return "operation is not implemented";
	case EFP_RESULT_ARGUMENT_NULL:
return "unexpected NULL argument to function";
	case EFP_RESULT_NOT_INITIALIZED:
return "efp was not properly initialized";
	case EFP_RESULT_FILE_NOT_FOUND:
return "file not found";
	case EFP_RESULT_SYNTAX_ERROR:
return "syntax error in potential data";
	case EFP_RESULT_UNKNOWN_FRAGMENT:
return "unknown fragment type";
	case EFP_RESULT_DUPLICATE_PARAMETERS:
return "fragment parameters contain fragments with the same name";
	case EFP_RESULT_CALLBACK_NOT_SET:
return "required callback function is not set";
	case EFP_RESULT_CALLBACK_FAILED:
return "callback function failed";
	case EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED:
return "overlap based damping requires exchange repulsion";
	case EFP_RESULT_GRADIENT_NOT_REQUESTED:
return "gradient computation was not requested";
	case EFP_RESULT_PARAMETERS_MISSING:
return "required EFP fragment parameters are missing";
	case EFP_RESULT_INDEX_OUT_OF_RANGE:
return "index is out of range";
	case EFP_RESULT_INVALID_ARRAY_SIZE:
return "invalid array size";
	case EFP_RESULT_UNSUPPORTED_SCREEN:
return "unsupported SCREEN group found in EFP data";
	case EFP_RESULT_INCONSISTENT_TERMS:
return "inconsistent energy terms selected";
	}
return "unknown result";
}
