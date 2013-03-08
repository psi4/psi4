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

#ifndef LIBEFP_PRIVATE_H
#define LIBEFP_PRIVATE_H

#include "efp.h"
#include "int.h"
#include "swf.h"
#include "terms.h"

#define EFP_EXPORT __attribute__((visibility("default")))
#define UNUSED __attribute__((unused))

#define ARRAY_SIZE(arr) (sizeof(arr)/sizeof(arr[0]))

struct frag {
	/* fragment name */
	char name[32];

	/* fragment center of mass */
	double x, y, z;

	/* rotation matrix representing orientation of a fragment */
	mat_t rotmat;

	/* pointer to the initial fragment state in library */
	const struct frag *lib;

	/* force on fragment center of mass */
	vec_t force;

	/* torque on fragment */
	vec_t torque;

	/* number of atoms in this fragment */
	int n_atoms;

	/* fragment atoms */
	struct efp_atom *atoms;

	/* distributed multipoles */
	struct multipole_pt {
		double x, y, z;
		double monopole;
		vec_t dipole;
		double quadrupole[6];
		double octupole[10];
	} *multipole_pts;

	/* number of distributed multipole points */
	int n_multipole_pts;

	/* electrostatic screening parameters */
	double *screen_params;

	/* ab initio electrostatic screening parameters */
	double *ai_screen_params;

	/* distributed polarizability points */
	struct polarizable_pt {
		double x, y, z;
		mat_t tensor;
		vec_t elec_field;
		vec_t elec_field_wf;
		vec_t induced_dipole;
		vec_t induced_dipole_new;
		vec_t induced_dipole_conj;
		vec_t induced_dipole_conj_new;
	} *polarizable_pts;

	/* number of distributed polarizability points */
	int n_polarizable_pts;

	/* dynamic polarizability points */
	struct dynamic_polarizable_pt {
		double x, y, z;
		double trace[12];
	} *dynamic_polarizable_pts;

	/* number of dynamic polarizability points */
	int n_dynamic_polarizable_pts;

	/* number of localized molecular orbitals */
	int n_lmo;

	/* localized molecular orbital centroids */
	vec_t *lmo_centroids;

	/* spin multiplicity */
	int multiplicity;

	/* number of exchange repulsion basis shells */
	int n_xr_shells;

	/* exchange repulsion basis shells */
	struct shell *xr_shells;

	/* upper triangle of fock matrix, size = n_lmo * (n_lmo + 1) / 2 */
	double *xr_fock_mat;

	/* exchange repulsion wavefunction size */
	int xr_wf_size;

	/* exchange repulsion wavefunction, size = n_lmo * xr_wf_size */
	double *xr_wf;

	/* rotational derivatives of MO coefficients */
	double *xr_wf_deriv[3];

	/* overlap integrals; used for overlap-based dispersion damping */
	double *overlap_int;

	/* derivatives of overlap integrals */
	six_t *overlap_int_deriv;
};

struct efp {
	/* number of fragments */
	int n_frag;

	/* array of fragments */
	struct frag *frags;

	/* number of fragments in the library */
	int n_lib;

	/* array with the library of fragment initial parameters */
	struct frag **lib;

	/* callback which computes electric field from electrons */
	efp_electron_density_field_fn get_electron_density_field;

	/* user data for get_electron_density_field */
	void *get_electron_density_field_user_data;

	/* user parameters for this EFP computation */
	struct efp_opts opts;

	/* gradient will also be computed if nonzero */
	int do_gradient;

	/* periodic simulation box size */
	vec_t box;

	/* stress tensor */
	mat_t stress;

	/* number of point charges */
	int n_ptc;

	struct point_charge {
		double x, y, z;
		double charge;
		vec_t grad;
	} *point_charges;

	/* EFP energy terms */
	struct efp_energy energy;
};

int efp_skip_frag_pair(struct efp *, int, int);
struct swf efp_make_swf(struct efp *, const struct frag *, const struct frag *);
const struct frag *efp_find_lib(struct efp *, const char *);
void efp_add_stress(const vec_t *, const vec_t *, mat_t *);
void efp_add_force(struct frag *, const vec_t *, const vec_t *, const vec_t *);
void efp_sub_force(struct frag *, const vec_t *, const vec_t *, const vec_t *);
void efp_move_pt(const vec_t *, const mat_t *, const vec_t *, vec_t *);
void efp_rotate_t2(const mat_t *, const double *, double *);
void efp_rotate_t3(const mat_t *, const double *, double *);
int efp_strcasecmp(const char *, const char *);
int efp_strncasecmp(const char *, const char *, size_t);

#endif /* LIBEFP_PRIVATE_H */
