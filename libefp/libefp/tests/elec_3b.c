#include "test_common.h"
#include "geometry_3.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_points = geometry,
	.ref_energy = 0.0023592829,
	.opts = {
		.terms = EFP_TERM_ELEC,
		.elec_damp = EFP_ELEC_DAMP_OVERLAP
	}
};

DEFINE_TEST(elec_3b)
