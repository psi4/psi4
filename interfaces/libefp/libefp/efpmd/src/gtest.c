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

void sim_gtest(struct state *state);

static void test_vec(char label, size_t idx, double tol, const double *agrad,
		const double *ngrad)
{
	bool match = true;

	for (size_t i = 0; i < 3; i++)
		if (fabs(agrad[i] - ngrad[i]) > tol)
			match = false;

	msg("A %c%04zu  ", label, idx);
	print_vec(agrad);
	msg("\n");
	msg("N %c%04zu  ", label, idx);
	print_vec(ngrad);
	msg(match ? "  MATCH\n" : "  DOES NOT MATCH\n");
}

static void test_cgrad(struct state *state, const double *cgrad)
{
	double tol = cfg_get_double(state->cfg, "gtest_tol");
	double dstep = cfg_get_double(state->cfg, "num_step_dist");

	size_t n_charges;
	check_fail(efp_get_point_charge_count(state->efp, &n_charges));

	double znuc[n_charges], xyz[3 * n_charges];
	check_fail(efp_get_point_charge_values(state->efp, znuc));
	check_fail(efp_get_point_charge_coordinates(state->efp, xyz));

	for (size_t i = 0; i < n_charges; i++) {
		double ngrad[3];

		for (size_t j = 0; j < 3; j++) {
			double e1, e2;
			double coord = xyz[3 * i + j];

			xyz[3 * i + j] = coord - dstep;
			check_fail(efp_set_point_charges(state->efp, n_charges, znuc, xyz));
			compute_energy(state, 0);
			e1 = state->energy;

			xyz[3 * i + j] = coord + dstep;
			check_fail(efp_set_point_charges(state->efp, n_charges, znuc, xyz));
			compute_energy(state, 0);
			e2 = state->energy;

			xyz[3 * i + j] = coord;
			ngrad[j] = (e2 - e1) / (2.0 * dstep);
		}

		test_vec('Q', i + 1, tol, cgrad + 3 * i, ngrad);
	}

	check_fail(efp_set_point_charges(state->efp, n_charges, znuc, xyz));
}

static void test_fgrad(struct state *state, const double *fgrad)
{
	double tol = cfg_get_double(state->cfg, "gtest_tol");
	double dstep = cfg_get_double(state->cfg, "num_step_dist");
	double astep = cfg_get_double(state->cfg, "num_step_angle");

	size_t n_frags;
	check_fail(efp_get_frag_count(state->efp, &n_frags));

	double xyzabc[6 * n_frags];
	check_fail(efp_get_coordinates(state->efp, xyzabc));

	for (size_t i = 0; i < n_frags; i++) {
		double deriv[3], ngrad[6];

		for (size_t j = 0; j < 6; j++) {
			double e1, e2;
			double coord = xyzabc[6 * i + j];
			double step = j < 3 ? dstep : astep;

			xyzabc[6 * i + j] = coord - step;
			check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, xyzabc));
			compute_energy(state, false);
			e1 = state->energy;

			xyzabc[6 * i + j] = coord + step;
			check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, xyzabc));
			compute_energy(state, false);
			e2 = state->energy;

			xyzabc[6 * i + j] = coord;
			ngrad[j] = (e2 - e1) / (2.0 * step);
		}

		test_vec('F', i + 1, tol, fgrad + 6 * i, ngrad);
		efp_torque_to_derivative(xyzabc + 6 * i + 3, fgrad + 6 * i + 3, deriv);
		test_vec('D', i + 1, tol, deriv, ngrad + 3);
	}

	check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, xyzabc));
}

static void test_grad(struct state *state)
{
	size_t n_frags, n_charges;
	check_fail(efp_get_frag_count(state->efp, &n_frags));
	check_fail(efp_get_point_charge_count(state->efp, &n_charges));

	double fgrad[6 * n_frags];
	memcpy(fgrad, state->grad, n_frags * 6 * sizeof(double));

	if (n_charges > 0) {
		double cgrad[3 * n_charges];
		check_fail(efp_get_point_charge_gradient(state->efp, cgrad));
		test_cgrad(state, cgrad);
	}

	test_fgrad(state, fgrad);
}

static void test_energy(struct state *state)
{
	double eref, tol;

	eref = cfg_get_double(state->cfg, "ref_energy");
	tol = cfg_get_double(state->cfg, "gtest_tol");

	msg("%30s %16.10lf\n", "REFERENCE ENERGY", eref);
	msg("%30s %16.10lf", "COMPUTED ENERGY", state->energy);
	msg(fabs(eref - state->energy) < tol ? "  MATCH\n" : "  DOES NOT MATCH\n");
}

void sim_gtest(struct state *state)
{
	msg("GRADIENT TEST JOB\n\n\n");

	print_geometry(state->efp);
	compute_energy(state, 1);
	print_energy(state);
	test_energy(state);

	msg("\n\n    COMPUTING NUMERICAL GRADIENT\n\n");
	test_grad(state);
	msg("\n");

	msg("GRADIENT TEST JOB COMPLETED SUCCESSFULLY\n");
}
