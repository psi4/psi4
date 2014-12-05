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

#include <math_util.h>

#include "clapack.h"
#include "common.h"

void sim_hess(struct state *state);

static void compute_gradient(struct state *state, size_t n_frags,
			     const double *xyzabc, double *grad)
{
	check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, xyzabc));
	compute_energy(state, true);
	memcpy(grad, state->grad, n_frags * 6 * sizeof(double));

	for (size_t i = 0; i < n_frags; i++) {
		const double *euler = xyzabc + 6 * i + 3;
		double *gradptr = grad + 6 * i + 3;

		efp_torque_to_derivative(euler, gradptr, gradptr);
	}
}

static void show_progress(size_t disp, size_t total, const char *dir)
{
	msg("COMPUTING DISPLACEMENT %4zu OF %zu (%s)\n", disp, total, dir);
	fflush(stdout);
}

static void compute_hessian(struct state *state, double *hess)
{
	size_t n_frags, n_coord;
	double *xyzabc, *grad_f, *grad_b;
	bool central = cfg_get_bool(state->cfg, "hess_central");

	check_fail(efp_get_frag_count(state->efp, &n_frags));
	n_coord = 6 * n_frags;

	xyzabc = xmalloc(n_coord * sizeof(double));
	grad_f = xmalloc(n_coord * sizeof(double));
	grad_b = xmalloc(n_coord * sizeof(double));

	check_fail(efp_get_coordinates(state->efp, xyzabc));

	if (!central) {
		memcpy(grad_b, state->grad, n_frags * 6 * sizeof(double));

		for (size_t i = 0; i < n_frags; i++) {
			const double *euler = xyzabc + 6 * i + 3;
			double *gradptr = grad_b + 6 * i + 3;

			efp_torque_to_derivative(euler, gradptr, gradptr);
		}
	}

	for (size_t i = 0; i < n_coord; i++) {
		double save = xyzabc[i];
		double step = i % 6 < 3 ? cfg_get_double(state->cfg, "num_step_dist") :
					  cfg_get_double(state->cfg, "num_step_angle");

		show_progress(i + 1, n_coord, "FORWARD");
		xyzabc[i] = save + step;
		compute_gradient(state, n_frags, xyzabc, grad_f);

		if (central) {
			show_progress(i + 1, n_coord, "BACKWARD");
			xyzabc[i] = save - step;
			compute_gradient(state, n_frags, xyzabc, grad_b);
		}

		double delta = central ? 2.0 * step : step;

		for (size_t j = 0; j < n_coord; j++)
			hess[i * n_coord + j] = (grad_f[j] - grad_b[j]) / delta;

		xyzabc[i] = save;
	}

	/* restore original coordinates */
	check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, xyzabc));

	/* reduce error by computing the average of H(i,j) and H(j,i) */
	for (size_t i = 0; i < n_coord; i++) {
		for (size_t j = i + 1; j < n_coord; j++) {
			double sum = hess[i * n_coord + j] + hess[j * n_coord + i];

			hess[i * n_coord + j] = 0.5 * sum;
			hess[j * n_coord + i] = hess[i * n_coord + j];
		}
	}

	free(xyzabc);
	free(grad_f);
	free(grad_b);

	msg("\n\n");
}

static void get_inertia_factor(const double *inertia, const mat_t *rotmat,
					mat_t *inertia_fact)
{
	double fact[3];

	for (size_t i = 0; i < 3; i++)
		fact[i] = inertia[i] < EPSILON ? 0.0 : 1.0 / sqrt(inertia[i]);

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			double sum = 0.0;

			for (size_t k = 0; k < 3; k++)
				sum += fact[k] * mat_get(rotmat, i, k) * mat_get(rotmat, j, k);

			mat_set(inertia_fact, i, j, sum);
		}
	}
}

static void get_weight_factor(struct efp *efp, double *mass_fact, mat_t *inertia_fact)
{
	size_t n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	double xyzabc[6 * n_frags];
	check_fail(efp_get_coordinates(efp, xyzabc));

	for (size_t i = 0; i < n_frags; i++) {
		double mass;
		check_fail(efp_get_frag_mass(efp, i, &mass));

		double inertia[3];
		check_fail(efp_get_frag_inertia(efp, i, inertia));

		double a = xyzabc[6 * i + 3];
		double b = xyzabc[6 * i + 4];
		double c = xyzabc[6 * i + 5];

		mat_t rotmat;
		euler_to_matrix(a, b, c, &rotmat);

		mass_fact[i] = 1.0 / sqrt(mass);
		get_inertia_factor(inertia, &rotmat, inertia_fact + i);
	}
}

static void w_tr_tr(double fact1, double fact2, size_t stride,
			const double *in, double *out)
{
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
			out[stride * i + j] = fact1 * fact2 * in[stride * i + j];
}

static void w_tr_rot(double fact1, const mat_t *fact2, size_t stride,
			const double *in, double *out)
{
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			double w = 0.0;

			for (size_t k = 0; k < 3; k++)
				w += mat_get(fact2, j, k) * in[stride * i + k];

			out[stride * i + j] = w * fact1;
		}
	}
}

static void w_rot_tr(const mat_t *fact1, double fact2, size_t stride,
			const double *in, double *out)
{
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			double w = 0.0;

			for (size_t k = 0; k < 3; k++)
				w += mat_get(fact1, i, k) * in[stride * k + j];

			out[stride * i + j] = w * fact2;
		}
	}
}

static void w_rot_rot(const mat_t *fact1, const mat_t *fact2, size_t stride,
			const double *in, double *out)
{
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			double w1 = 0.0;

			for (size_t ii = 0; ii < 3; ii++) {
				double w2 = 0.0;

				for (size_t jj = 0; jj < 3; jj++)
					w2 += mat_get(fact2, j, jj) * in[stride * ii + jj];

				w1 += w2 * mat_get(fact1, i, ii);
			}

			out[i * stride + j] = w1;
		}
	}
}

static void mass_weight_hessian(struct efp *efp, const double *in, double *out)
{
	size_t n_frags, n_coord;

	check_fail(efp_get_frag_count(efp, &n_frags));
	n_coord = 6 * n_frags;

	double mass_fact[n_frags];
	mat_t inertia_fact[n_frags];

	get_weight_factor(efp, mass_fact, inertia_fact);

	for (size_t i = 0; i < n_frags; i++) {
		for (size_t j = 0; j < n_frags; j++) {
			size_t offset = 6 * n_coord * i + 6 * j;

			w_tr_tr(mass_fact[i], mass_fact[j], n_coord,
					in + offset,
					out + offset);

			w_tr_rot(mass_fact[i], inertia_fact + j, n_coord,
					in + offset + 3,
					out + offset + 3);

			w_rot_tr(inertia_fact + i, mass_fact[j], n_coord,
					in + offset + 3 * n_coord,
					out + offset + 3 * n_coord);

			w_rot_rot(inertia_fact + i, inertia_fact + j, n_coord,
					in + offset + 3 * n_coord + 3,
					out + offset + 3 * n_coord + 3);
		}
	}
}

static void print_mode(size_t mode, double eigen)
{
	/* preserve sign for imaginary frequencies */
	eigen = copysign(sqrt(fabs(eigen)), eigen);

	/* convert to cm-1 */
	eigen = FINE_CONST / 2.0 / PI / BOHR_RADIUS / sqrt(AMU_TO_AU) * 1.0e8 * eigen;

	msg("    MODE %4zu    FREQUENCY %10.3lf cm-1\n\n", mode, eigen);
}

void sim_hess(struct state *state)
{
	msg("HESSIAN JOB\n\n\n");

	print_geometry(state->efp);
	compute_energy(state, true);
	print_energy(state);
	print_gradient(state);

	size_t n_frags, n_coord;
	double *hess, *mass_hess, *eigen;

	check_fail(efp_get_frag_count(state->efp, &n_frags));
	n_coord = 6 * n_frags;

	hess = xmalloc(n_coord * n_coord * sizeof(double));
	compute_hessian(state, hess);

	msg("    HESSIAN MATRIX\n\n");
	print_matrix(n_coord, n_coord, hess);

	mass_hess = xmalloc(n_coord * n_coord * sizeof(double));
	mass_weight_hessian(state->efp, hess, mass_hess);

	msg("    MASS-WEIGHTED HESSIAN MATRIX\n\n");
	print_matrix(n_coord, n_coord, mass_hess);

	msg("    NORMAL MODE ANALYSIS\n\n");

	eigen = xmalloc(n_coord * sizeof(double));

	if (efp_dsyev('V', 'U', (int)n_coord, mass_hess, (int)n_coord, eigen))
		error("unable to diagonalize mass-weighted hessian matrix");

	for (size_t i = 0; i < n_coord; i++) {
		print_mode(i + 1, eigen[i]);
		print_vector(n_coord, mass_hess + i * n_coord);
	}

	free(hess);
	free(mass_hess);
	free(eigen);

	msg("HESSIAN JOB COMPLETED SUCCESSFULLY\n");
}
