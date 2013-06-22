#include "test_common.h"
#include "geometry_3.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_points = geometry,
	.ref_energy = -0.0173897265,
	.opts = {
		.terms = EFP_TERM_DISP,
		.disp_damp = EFP_DISP_DAMP_TT
	}
};

DEFINE_TEST(disp_3a)
