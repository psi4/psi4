#include "test_common.h"
#include "geometry_4.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = -0.0092400662,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR,
		.elec_damp = EFP_ELEC_DAMP_OVERLAP,
		.disp_damp = EFP_DISP_DAMP_OVERLAP
	}
};

DEFINE_TEST(total_4b)
