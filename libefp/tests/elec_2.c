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
	-0.000161457,   -0.000684540,   -0.001599534,
	 0.000900919,   -0.000053204,   -0.002052563,
	-0.000090312,    0.000134507,   -0.000385729,
	-0.000052205,   -0.000612205,   -0.000080863,
	 0.000443897,    0.000339330,    0.000687844,
	 0.001232655,    0.000092766,   -0.000462341,
	-0.000079313,    0.000322287,    0.000933484,
	-0.001671896,    0.000461100,    0.000512841,
	-0.000112816,   -0.000111583,    0.000363934,
	 0.000136646,   -0.000779130,   -0.001022605
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.geometry_xyzabc = xyzabc,
	.ref_energy = 0.0013721463, /* from GAMESS */
	.energy_accuracy = 6,
	.ref_gradient = ref_gradient,
	.gradient_accuracy = 6,
	.test_numerical_gradient = 1,
	.opts = {
		.terms = EFP_TERM_ELEC
	}
};

DEFINE_TEST(test_data)
