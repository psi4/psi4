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

#include "common.h"
#include "opt.h"

void sim_opt(struct state *state);

static double compute_efp(size_t n, const double *x, double *gx, void *data)
{
	size_t n_frags, n_charge;
	struct state *state = (struct state *)data;

	check_fail(efp_get_frag_count(state->efp, &n_frags));
	check_fail(efp_get_point_charge_count(state->efp, &n_charge));

	assert(n == (6 * n_frags + 3 * n_charge));

	check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, x));
	check_fail(efp_set_point_charge_coordinates(state->efp, x + 6 * n_frags));

	compute_energy(state, true);
	memcpy(gx, state->grad, (6 * n_frags + 3 * n_charge) * sizeof(double));

	for (size_t i = 0; i < n_frags; i++) {
		const double *euler = x + 6 * i + 3;
		double *gradptr = gx + 6 * i + 3;

		efp_torque_to_derivative(euler, gradptr, gradptr);
	}

	return (state->energy);
}

static void print_restart(struct efp *efp)
{
	size_t n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	double coord[6 * n_frags];
	check_fail(efp_get_coordinates(efp, coord));

	msg("    RESTART DATA\n\n");

	for (size_t i = 0; i < n_frags; i++) {
		char name[64];
		check_fail(efp_get_frag_name(efp, i, sizeof(name), name));

		coord[6 * i + 0] *= BOHR_RADIUS;
		coord[6 * i + 1] *= BOHR_RADIUS;
		coord[6 * i + 2] *= BOHR_RADIUS;

		print_fragment(name, coord + 6 * i, NULL);
	}

	size_t n_charges;
	check_fail(efp_get_point_charge_count(efp, &n_charges));

	if (n_charges > 0) {
		double q[n_charges];
		check_fail(efp_get_point_charge_values(efp, q));

		double xyz[3 * n_charges];
		check_fail(efp_get_point_charge_coordinates(efp, xyz));

		for (size_t i = 0; i < n_charges; i++) {
			double x = xyz[3 * i + 0] * BOHR_RADIUS;
			double y = xyz[3 * i + 1] * BOHR_RADIUS;
			double z = xyz[3 * i + 2] * BOHR_RADIUS;

			print_charge(q[i], x, y, z);
		}
	}

	msg("\n");
}

static int check_conv(double rms_grad, double max_grad, double opt_tol)
{
	return max_grad < opt_tol && rms_grad < opt_tol / 3.0;
}

static void get_grad_info(size_t n_coord, const double *grad, double *rms_grad_out,
				double *max_grad_out)
{
	double rms_grad = 0.0, max_grad = 0.0;

	for (size_t i = 0; i < n_coord; i++) {
		rms_grad += grad[i] * grad[i];

		if (fabs(grad[i]) > max_grad)
			max_grad = fabs(grad[i]);
	}

	rms_grad = sqrt(rms_grad / n_coord);

	*rms_grad_out = rms_grad;
	*max_grad_out = max_grad;
}

static void print_status(struct state *state, double e_diff, double rms_grad, double max_grad)
{
	print_geometry(state->efp);
	print_restart(state->efp);
	print_energy(state);

	msg("%30s %16.10lf\n", "ENERGY CHANGE", e_diff);
	msg("%30s %16.10lf\n", "RMS GRADIENT", rms_grad);
	msg("%30s %16.10lf\n", "MAXIMUM GRADIENT", max_grad);
	msg("\n\n");

	fflush(stdout);
}

void sim_opt(struct state *state)
{
	msg("ENERGY MINIMIZATION JOB\n\n\n");

	size_t n_frags, n_charge, n_coord;
	double rms_grad, max_grad;

	check_fail(efp_get_frag_count(state->efp, &n_frags));
	check_fail(efp_get_point_charge_count(state->efp, &n_charge));

	n_coord = 6 * n_frags + 3 * n_charge;

	struct opt_state *opt_state = opt_create(n_coord);
	if (!opt_state)
		error("unable to create an optimizer");

	opt_set_func(opt_state, compute_efp);
	opt_set_user_data(opt_state, state);

	double coord[n_coord], grad[n_coord];
	check_fail(efp_get_coordinates(state->efp, coord));
	check_fail(efp_get_point_charge_coordinates(state->efp, coord + 6 * n_frags));

	if (opt_init(opt_state, n_coord, coord))
		error("unable to initialize an optimizer");

	double e_old = opt_get_fx(opt_state);
	opt_get_gx(opt_state, n_coord, grad);
	get_grad_info(n_coord, grad, &rms_grad, &max_grad);

	msg("    INITIAL STATE\n\n");
	print_status(state, 0.0, rms_grad, max_grad);

	for (int step = 1; step <= cfg_get_int(state->cfg, "max_steps"); step++) {
		if (opt_step(opt_state))
			error("unable to make an optimization step");

		double e_new = opt_get_fx(opt_state);
		opt_get_gx(opt_state, n_coord, grad);
		get_grad_info(n_coord, grad, &rms_grad, &max_grad);

		if (check_conv(rms_grad, max_grad, cfg_get_double(state->cfg, "opt_tol"))) {
			msg("    FINAL STATE\n\n");
			print_status(state, e_new - e_old, rms_grad, max_grad);
			msg("OPTIMIZATION CONVERGED\n");
			break;
		}

		msg("    STATE AFTER %d STEPS\n\n", step);
		print_status(state, e_new - e_old, rms_grad, max_grad);

		e_old = e_new;
	}

	opt_shutdown(opt_state);

	msg("ENERGY MINIMIZATION JOB COMPLETED SUCCESSFULLY\n");
}
