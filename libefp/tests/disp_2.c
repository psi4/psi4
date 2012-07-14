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

static const double ref_gradient[] = { /* from Q-Chem 4.0 */
    -7.8926221726694876e-05,  0.0001952345800678e+00,  0.0001609972134529e+00,
     1.4638347286741127e-06,  1.3542830761268802e-06, -3.4820011499185243e-06,
     0.0004075855774725e+00, -0.0002238688386432e+00, -4.0220174676333711e-05,
     2.2661761976782832e-06, -7.1426825847379459e-06,  5.2255276600075696e-05,
     9.2974601898653485e-05, -1.1388610209512631e-05,  0.0001298654246440e+00,
     2.0112727120278625e-06,  1.4907395484342012e-06, -1.2395432287061943e-06,
     5.9230965170915758e-05,  5.2711048865655091e-05, -0.0002240681861481e+00,
     2.4255180334939735e-06, -2.0547661655446631e-06, -7.4874526498860911e-07,
    -0.0004808649228154e+00, -1.2688180080797921e-05, -2.6574277272524212e-05,
     9.0883272966491658e-06,  2.8159597899192419e-05, -1.3091016031986133e-05
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.geometry_xyzabc = xyzabc,
	.ref_energy = -0.0014688094, /* from GAMESS */
	.ref_gradient = ref_gradient,
	.test_numerical_gradient = 1,
	.opts = {
		.terms = EFP_TERM_DISP,
		.disp_damp = EFP_DISP_DAMP_TT
	}
};

DEFINE_TEST(test_data)
