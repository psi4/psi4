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

#ifndef LIBEFP_TEST_COMMON_H
#define LIBEFP_TEST_COMMON_H

#include <stdio.h>

#include <check.h>
#include <efp.h>

#include "test_list.h"

#define BOHR_RADIUS 0.52917721092
#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))
#define BOHR(x) ((x) / BOHR_RADIUS)
#define ANGSTROM(x) ((x) * BOHR_RADIUS)

struct test_data {
	/** Paths to EFP data files. */
	const char **files;

	/** Fragment names. */
	const char **names;

	/**
	 * If not NULL geometry will be set up using EFP_COORD_TYPE_XYZABC. */
	const double *geometry_xyzabc;

	/**
	 * If not NULL geometry will be set up using EFP_COORD_TYPE_POINTS. */
	const double *geometry_points;

	/** Periodic box size. */
	double box[3];

	/** Number of point charges. */
	int n_ptc;

	/** Point charges. */
	const double *ptc_charges;

	/** Coordinates of point charges. */
	const double *ptc_xyz;

	/** Reference energy value. */
	double ref_energy;

	/** Simulation settings. */
	struct efp_opts opts;

	/** Callback which computes electric field from electrons. */
	efp_electron_density_field_fn electron_density_field_fn;

	/** User data for electron_density_field_fn. */
	void *electron_density_field_user_data;
};

void run_test(const struct test_data *);

#define DEFINE_TEST(name)                                                    \
START_TEST(test)                                                             \
{                                                                            \
	printf("Running test %s\n", #name);                                  \
	run_test(&test_data);                                                \
}                                                                            \
END_TEST                                                                     \
TCase *tcase_ ## name(void)                                                  \
{                                                                            \
	TCase *tc = tcase_create(#name);                                     \
	tcase_add_test(tc, test);                                            \
	tcase_set_timeout(tc, 100);                                          \
	return tc;                                                           \
}

#endif /* LIBEFP_TEST_COMMON_H */
