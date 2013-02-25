#include "test_common.h"
#include "geometry_2.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = 0.0013685212,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL
	}
};

DEFINE_TEST(pol_2)
