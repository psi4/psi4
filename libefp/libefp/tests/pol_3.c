#include "test_common.h"
#include "geometry_3.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_points = geometry,
	.ref_energy = -0.0066095992,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL
	}
};

DEFINE_TEST(pol_3)
