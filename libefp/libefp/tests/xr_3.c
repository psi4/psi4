#include "test_common.h"
#include "geometry_3.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_points = geometry,
	.ref_energy = 0.0301402098,
	.opts = {
		.terms = EFP_TERM_XR
	}
};

DEFINE_TEST(xr_3)
