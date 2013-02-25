#include "test_common.h"
#include "geometry_2.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = 0.0015865516,
	.opts = {
		.terms = EFP_TERM_ELEC
	}
};

DEFINE_TEST(elec_2a)
