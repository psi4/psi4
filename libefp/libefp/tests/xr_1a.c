#include "test_common.h"
#include "geometry_1.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = 0.0000134716,
	.opts = {
		.terms = EFP_TERM_XR
	}
};

DEFINE_TEST(xr_1a)
