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

#include "common.h"
#include "optimizer.h"

void sim_opt(struct efp *, const struct config *);

static double compute_efp(int n, const double *x, double *gx, void *data)
{
	struct efp *efp = (struct efp *)data;

	int n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	assert(n == 6 * n_frags);

	check_fail(efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, x));
	check_fail(efp_compute(efp, 1));

	struct efp_energy energy;
	check_fail(efp_get_energy(efp, &energy));

	check_fail(efp_get_gradient(efp, n_frags, gx));

	for (int i = 0; i < n_frags; i++) {
		const double *euler = x + 6 * i + 3;
		double *gradptr = gx + 6 * i + 3;

		efp_torque_to_derivative(euler, gradptr, gradptr);
	}

	return energy.total;
}

static void print_restart(struct efp *efp)
{
	int n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	double coord[6 * n_frags];
	check_fail(efp_get_coordinates(efp, n_frags, coord));

	printf("    RESTART DATA\n\n");

	for (int i = 0; i < n_frags; i++) {
		char name[64];
		check_fail(efp_get_frag_name(efp, i, sizeof(name), name));

		coord[6 * i + 0] *= BOHR_RADIUS;
		coord[6 * i + 1] *= BOHR_RADIUS;
		coord[6 * i + 2] *= BOHR_RADIUS;

		print_fragment(name, coord + 6 * i, NULL);
	}

	printf("\n");
}

static int check_conv(double rms_grad, double max_grad, double opt_tol)
{
	return max_grad < opt_tol && rms_grad < opt_tol / 3.0;
}

static void get_grad_info(int n_coord, const double *grad, double *rms_grad_out,
				double *max_grad_out)
{
	double rms_grad = 0.0, max_grad = 0.0;

	for (int i = 0; i < n_coord; i++) {
		rms_grad += grad[i] * grad[i];

		if (fabs(grad[i]) > max_grad)
			max_grad = fabs(grad[i]);
	}

	rms_grad = sqrt(rms_grad / n_coord);

	*rms_grad_out = rms_grad;
	*max_grad_out = max_grad;
}

static void print_status(struct efp *efp, double e_diff, double rms_grad,
				double max_grad)
{
	print_geometry(efp);
	print_restart(efp);
	print_energy(efp);

	printf("%30s %16.10lf\n", "ENERGY CHANGE", e_diff);
	printf("%30s %16.10lf\n", "RMS GRADIENT", rms_grad);
	printf("%30s %16.10lf\n", "MAXIMUM GRADIENT", max_grad);
	printf("\n\n");

	fflush(stdout);
}

void sim_opt(struct efp *efp, const struct config *config)
{
	printf("ENERGY MINIMIZATION JOB\n\n\n");

	int n_coord = 6 * config->n_frags;
	double rms_grad, max_grad;

	struct opt_state *state = opt_create(n_coord);
	if (!state)
		error("UNABLE TO CREATE AN OPTIMIZER");

	opt_set_func(state, compute_efp);
	opt_set_user_data(state, efp);

	double coord[n_coord], grad[n_coord];
	check_fail(efp_get_coordinates(efp, config->n_frags, coord));

	if (opt_init(state, n_coord, coord))
		error("UNABLE TO INITIALIZE AN OPTIMIZER");

	double e_old = opt_get_fx(state);
	opt_get_gx(state, n_coord, grad);
	get_grad_info(n_coord, grad, &rms_grad, &max_grad);

	printf("    INITIAL STATE\n\n");
	print_status(efp, 0.0, rms_grad, max_grad);

	for (int step = 1; step <= config->max_steps; step++) {
		if (opt_step(state))
			error("UNABLE TO MAKE AN OPTIMIZATION STEP");

		double e_new = opt_get_fx(state);
		opt_get_gx(state, n_coord, grad);
		get_grad_info(n_coord, grad, &rms_grad, &max_grad);

		if (check_conv(rms_grad, max_grad, config->opt_tol)) {
			printf("    FINAL STATE\n\n");
			print_status(efp, e_new - e_old, rms_grad, max_grad);
			printf("OPTIMIZATION CONVERGED\n");
			break;
		}

		printf("    STATE AFTER %d STEPS\n\n", step);
		print_status(efp, e_new - e_old, rms_grad, max_grad);

		e_old = e_new;
	}

	opt_shutdown(state);

	printf("ENERGY MINIMIZATION JOB COMPLETED SUCCESSFULLY\n");
}
