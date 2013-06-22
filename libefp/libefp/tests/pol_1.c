#include "test_common.h"
#include "geometry_1.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = 0.0002777238,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL
	}
};

DEFINE_TEST(pol_1)
