#include "test_common.h"
#include "geometry_3.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_points = geometry,
	.ref_energy = -0.0220107872,
	.opts = {
		.terms = EFP_TERM_DISP,
		.disp_damp = EFP_DISP_DAMP_OVERLAP
	}
};

DEFINE_TEST(disp_3b)
