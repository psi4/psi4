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

#include <math.h>
#include <stdlib.h>

#include "test_common.h"

#define FAIL_TOL 5.0e-6
#define NUM_GRAD_DELTA 0.001

static void lib_error(const char *title, enum efp_result res)
{
	fail("%s:\n    %s\n", title, efp_result_to_string(res));
}

static void test_qm_numerical_grad(struct efp *efp, const double *grad)
{
	int n_ptc;
	enum efp_result res;

	if ((res = efp_get_point_charge_count(efp, &n_ptc)))
		lib_error("efp_get_point_charge_count", res);

	double znuc[n_ptc], xyz[3 * n_ptc];

	if ((res = efp_get_point_charge_coordinates(efp, xyz)))
		lib_error("efp_get_point_charge_coordinates", res);

	if ((res = efp_get_point_charge_values(efp, znuc)))
		lib_error("efp_get_point_charge_values", res);

	for (int i = 0; i < n_ptc; i++) {
		double grad_num[3];

		for (int j = 0; j < 3; j++) {
			double coord = xyz[3 * i + j];

			struct efp_energy e1;
			xyz[3 * i + j] = coord - NUM_GRAD_DELTA;

			if ((res = efp_set_point_charges(efp, n_ptc, znuc, xyz)))
				lib_error("efp_set_point_charges", res);
			if ((res = efp_compute(efp, 0)))
				lib_error("efp_compute", res);
			if ((res = efp_get_energy(efp, &e1)))
				lib_error("efp_get_energy", res);

			struct efp_energy e2;
			xyz[3 * i + j] = coord + NUM_GRAD_DELTA;

			if ((res = efp_set_point_charges(efp, n_ptc, znuc, xyz)))
				lib_error("efp_set_point_charges", res);
			if ((res = efp_compute(efp, 0)))
				lib_error("efp_compute", res);
			if ((res = efp_get_energy(efp, &e2)))
				lib_error("efp_get_energy", res);

			xyz[3 * i + j] = coord;
			grad_num[j] = (e2.total - e1.total) / (2.0 * NUM_GRAD_DELTA);
		}

		for (int j = 0; j < 3; j++)
			fail_unless(fabs(grad_num[j] - grad[3 * i + j]) < FAIL_TOL);
	}

	if ((res = efp_set_point_charges(efp, n_ptc, znuc, xyz)))
		lib_error("efp_set_point_charges", res);
}

static void test_frag_numerical_grad(struct efp *efp, double *xyzabc, const double *grad)
{
	int n_frag;
	enum efp_result res;

	if ((res = efp_get_frag_count(efp, &n_frag)))
		lib_error("efp_get_frag_count", res);

	for (int i = 0; i < n_frag; i++) {
		double grad_num[6];

		for (int j = 0; j < 6; j++) {
			double coord = xyzabc[6 * i + j];

			struct efp_energy e1;
			xyzabc[6 * i + j] = coord - NUM_GRAD_DELTA;

			if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
				lib_error("efp_set_coordinates", res);
			if ((res = efp_compute(efp, 0)))
				lib_error("efp_compute", res);
			if ((res = efp_get_energy(efp, &e1)))
				lib_error("efp_get_energy", res);

			struct efp_energy e2;
			xyzabc[6 * i + j] = coord + NUM_GRAD_DELTA;

			if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
				lib_error("efp_set_coordinates", res);
			if ((res = efp_compute(efp, 0)))
				lib_error("efp_compute", res);
			if ((res = efp_get_energy(efp, &e2)))
				lib_error("efp_get_energy", res);

			xyzabc[6 * i + j] = coord;
			grad_num[j] = (e2.total - e1.total) / (2.0 * NUM_GRAD_DELTA);
		}

		for (int j = 0; j < 6; j++)
			fail_unless(fabs(grad_num[j] - grad[6 * i + j]) < FAIL_TOL);
	}

	if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
		lib_error("efp_set_coordinates", res);
}

void run_test(const struct test_data *test_data)
{
	enum efp_result res;
	struct efp *efp = efp_create();

	if (!efp)
		fail("efp_create failed\n");

	if ((res = efp_set_opts(efp, &test_data->opts)))
		lib_error("efp_set_opts", res);

	for (int i = 0; test_data->files[i]; i++)
		if ((res = efp_add_potential(efp, test_data->files[i])))
			lib_error("efp_add_potential", res);

	for (int i = 0; test_data->names[i]; i++)
		if ((res = efp_add_fragment(efp, test_data->names[i])))
			lib_error("efp_add_fragment", res);

	if ((res = efp_set_electron_density_field_fn(efp,
			test_data->electron_density_field_fn)))
		lib_error("efp_set_electron_density_field_fn", res);

	if ((res = efp_set_electron_density_field_user_data(efp,
			test_data->electron_density_field_user_data)))
		lib_error("efp_set_electron_density_field_user_data", res);

	if (test_data->geometry_xyzabc) {
		if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC,
						test_data->geometry_xyzabc)))
			lib_error("efp_set_coordinates", res);
	}
	else if (test_data->geometry_points) {
		if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_POINTS,
						test_data->geometry_points)))
			lib_error("efp_set_coordinates", res);
	}
	else {
		fail("geometry is not set");
	}

	int do_qm = test_data->ptc_charges && test_data->ptc_xyz;

	if (do_qm) {
		if ((res = efp_set_point_charges(efp, test_data->n_ptc,
					test_data->ptc_charges, test_data->ptc_xyz)))
			lib_error("efp_set_point_charges", res);
	}

	if (test_data->opts.enable_pbc) {
		double bx = test_data->box[0];
		double by = test_data->box[1];
		double bz = test_data->box[2];

		if ((res = efp_set_periodic_box(efp, bx, by, bz)))
			lib_error("efp_set_periodic_box", res);
	}

	/* Begin imaginary ab initio SCF */
	double scf_energy;

	if ((res = efp_get_wavefunction_dependent_energy(efp, &scf_energy)))
		lib_error("efp_get_wavefunction_dependent_energy", res);
	/* End imaginary ab initio SCF */

	if ((res = efp_compute(efp, 1)))
		lib_error("efp_compute", res);

	struct efp_energy energy;

	if ((res = efp_get_energy(efp, &energy)))
		lib_error("efp_get_energy", res);

	fail_unless(fabs(energy.total - test_data->ref_energy) < FAIL_TOL);

	int n_frag;
	if ((res = efp_get_frag_count(efp, &n_frag)))
		lib_error("efp_get_frag_count", res);

	double frag_grad[6 * n_frag];
	if ((res = efp_get_gradient(efp, n_frag, frag_grad)))
		lib_error("efp_get_gradient", res);

	double xyzabc[6 * n_frag];
	if ((res = efp_get_coordinates(efp, n_frag, xyzabc)))
		lib_error("efp_get_coordinates", res);

	for (int i = 0; i < n_frag; i++) {
		const double *euler = xyzabc + 6 * i + 3;
		double *gradptr = frag_grad + 6 * i + 3;

		efp_torque_to_derivative(euler, gradptr, gradptr);
	}

	if (do_qm) {
		int n_ptc;
		if ((res = efp_get_point_charge_count(efp, &n_ptc)))
			lib_error("efp_get_point_charge_count", res);

		double qm_grad[3 * n_ptc];
		if ((res = efp_get_point_charge_gradient(efp, qm_grad)))
			lib_error("efp_get_point_charge_gradient", res);

		test_qm_numerical_grad(efp, qm_grad);
	}

	test_frag_numerical_grad(efp, xyzabc, frag_grad);

	efp_shutdown(efp);
}

static void add_tcases(Suite *s)
{
	suite_add_tcase(s, tcase_disp_1a());
	suite_add_tcase(s, tcase_disp_1b());
	suite_add_tcase(s, tcase_disp_1c());
	suite_add_tcase(s, tcase_disp_1d());
	suite_add_tcase(s, tcase_disp_2a());
	suite_add_tcase(s, tcase_disp_2b());
	suite_add_tcase(s, tcase_disp_3a());
	suite_add_tcase(s, tcase_disp_3b());
	suite_add_tcase(s, tcase_elec_1a());
	suite_add_tcase(s, tcase_elec_1b());
	suite_add_tcase(s, tcase_elec_1c());
	suite_add_tcase(s, tcase_elec_2a());
	suite_add_tcase(s, tcase_elec_2b());
	suite_add_tcase(s, tcase_elec_3a());
	suite_add_tcase(s, tcase_elec_3b());
	suite_add_tcase(s, tcase_pol_1());
	suite_add_tcase(s, tcase_pol_2());
	suite_add_tcase(s, tcase_pol_3());
	suite_add_tcase(s, tcase_qm_1());
	suite_add_tcase(s, tcase_qm_2());
	suite_add_tcase(s, tcase_total_1a());
	suite_add_tcase(s, tcase_total_2a());
	suite_add_tcase(s, tcase_total_3a());
	suite_add_tcase(s, tcase_total_4a());
	suite_add_tcase(s, tcase_total_4b());
	suite_add_tcase(s, tcase_total_4c());
	suite_add_tcase(s, tcase_total_5a());
	suite_add_tcase(s, tcase_total_5b());
	suite_add_tcase(s, tcase_total_5c());
	suite_add_tcase(s, tcase_total_6a());
	suite_add_tcase(s, tcase_total_6b());
	suite_add_tcase(s, tcase_total_6c());
	suite_add_tcase(s, tcase_xr_1a());
	suite_add_tcase(s, tcase_xr_1b());
	suite_add_tcase(s, tcase_xr_2());
	suite_add_tcase(s, tcase_xr_3());
}

int main(void)
{
	Suite *s = suite_create("libefp test suite");

	add_tcases(s);

	SRunner *sr = srunner_create(s);

	srunner_set_fork_status(sr, CK_FORK);
	srunner_run_all(sr, CK_NORMAL);

	int n_failed = srunner_ntests_failed(sr);

	srunner_free(sr);

	return n_failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
