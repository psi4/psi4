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

#ifndef LIBEFP_EFP_PRIVATE_H
#define LIBEFP_EFP_PRIVATE_H

#include <assert.h>

#include "efp.h"
#include "int.h"
#include "math_util.h"
#include "terms.h"

#define EFP_EXPORT __attribute__((visibility("default")))
#define ARRAY_SIZE(arr) (sizeof(arr)/sizeof(arr[0]))
#define EFP_INIT_MAGIC 0xEF2012AD

struct frag {
	/* fragment name */
	char *name;

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
};

struct efp {
	/* number of fragments */
	int n_frag;

	/* array of fragments */
	struct frag *frags;

	/* number of fragments in the library */
	int n_lib;

	/* array with the library of fragment initial parameters */
	struct frag *lib;

	/* contributions to dispersion damping from overlap integrals */
	double *disp_damp_overlap;

	/* fragment offsets in disp_damp_overlap array */
	int *disp_damp_overlap_offset;

	/* callbacks */
	struct efp_callbacks callbacks;

	/* user parameters for this EFP computation */
	struct efp_opts opts;

	/* gradient will also be computed if nonzero */
	int do_gradient;

	/* information about ab initio region */
	struct {
		int n_atoms;
		double *znuc;
		double *xyz;
		double *grad;
	} qm;

	/* EFP energy terms */
	struct efp_energy energy;

	/* initialization check */
	unsigned magic;
};

static inline void
add_force_torque(struct frag *fr_i, struct frag *fr_j,
		 const vec_t *pt_i, const vec_t *pt_j,
		 const vec_t *force)
{
	vec_t dr_i = vec_sub(VEC(pt_i->x), VEC(fr_i->x));
	vec_t dr_j = vec_sub(VEC(pt_j->x), VEC(fr_j->x));

	vec_t torque_i = vec_cross(&dr_i, force);
	vec_t torque_j = vec_cross(&dr_j, force);

	#pragma omp atomic
	fr_i->force.x += force->x;
	#pragma omp atomic
	fr_i->force.y += force->y;
	#pragma omp atomic
	fr_i->force.z += force->z;

	#pragma omp atomic
	fr_i->torque.x += torque_i.x;
	#pragma omp atomic
	fr_i->torque.y += torque_i.y;
	#pragma omp atomic
	fr_i->torque.z += torque_i.z;

	#pragma omp atomic
	fr_j->force.x -= force->x;
	#pragma omp atomic
	fr_j->force.y -= force->y;
	#pragma omp atomic
	fr_j->force.z -= force->z;

	#pragma omp atomic
	fr_j->torque.x -= torque_j.x;
	#pragma omp atomic
	fr_j->torque.y -= torque_j.y;
	#pragma omp atomic
	fr_j->torque.z -= torque_j.z;
}

static inline void
add_force_torque_2(struct frag *fr_i, struct frag *fr_j,
		   const vec_t *pt_i, const vec_t *pt_j,
		   const vec_t *force, const vec_t *add_i, const vec_t *add_j)
{
	add_force_torque(fr_i, fr_j, pt_i, pt_j, force);

	#pragma omp atomic
	fr_i->torque.x += add_i->x;
	#pragma omp atomic
	fr_i->torque.y += add_i->y;
	#pragma omp atomic
	fr_i->torque.z += add_i->z;

	#pragma omp atomic
	fr_j->torque.x -= add_j->x;
	#pragma omp atomic
	fr_j->torque.y -= add_j->y;
	#pragma omp atomic
	fr_j->torque.z -= add_j->z;
}

#endif /* LIBEFP_EFP_PRIVATE_H */
