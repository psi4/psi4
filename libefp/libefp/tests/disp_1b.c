#include "test_common.h"
#include "geometry_1.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = -0.0001007275,
	.opts = {
		.terms = EFP_TERM_DISP,
		.disp_damp = EFP_DISP_DAMP_OVERLAP
	}
};

DEFINE_TEST(disp_1b)
