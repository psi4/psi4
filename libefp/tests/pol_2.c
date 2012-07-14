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

#include "test_common.h"
#include "geometry_2.h"

static const double ref_gradient[] = { /* from GAMESS */
	-0.000168166,   -0.000676699,   -0.001532420,
	 0.000897361,   -0.000032039,   -0.002046575,
	-0.000060465,    0.000108415,   -0.000381347,
	-0.000081349,   -0.000622769,   -0.000032799,
	 0.000454457,    0.000345412,    0.000717751,
	 0.001241851,    0.000122818,   -0.000498105,
	-0.000077235,    0.000327348,    0.000855399,
	-0.001629347,    0.000422712,    0.000506788,
	-0.000148591,   -0.000104477,    0.000340617,
	 0.000154698,   -0.000709890,   -0.001015295
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.geometry_xyzabc = xyzabc,
		/* elec + pol - from GAMESS */
	.ref_energy = 0.0013721463 + -0.0001902044,
	.energy_accuracy = 6,
	.ref_gradient = ref_gradient,
	.gradient_accuracy = 5,
	.test_numerical_gradient = 1,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL
	}
};

DEFINE_TEST(test_data)
