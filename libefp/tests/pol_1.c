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
#include "geometry_1.h"

static const double ref_gradient[] = { /* from GAMESS */
	 0.000096311,   -0.000039096,   -0.000046666,
	 0.000060285,    0.000530526,   -0.000252051,
	-0.000096311,    0.000039096,    0.000046666,
	-0.000060285,   -0.000089595,   -0.000117352
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.geometry_xyzabc = xyzabc,
		/* elec + pol - from GAMESS */
	.ref_energy = 0.0002554254 + -0.0000104028,
	.ref_gradient = ref_gradient,
	.gradient_accuracy = 6,
	.test_numerical_gradient = 1,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL
	}
};

DEFINE_TEST(test_data)
