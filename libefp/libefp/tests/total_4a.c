#include "test_common.h"
#include "geometry_4.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = -0.0095597483,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR,
		.elec_damp = EFP_ELEC_DAMP_SCREEN,
		.disp_damp = EFP_DISP_DAMP_TT
	}
};

DEFINE_TEST(total_4a)
