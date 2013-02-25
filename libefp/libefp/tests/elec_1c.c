#include "test_common.h"
#include "geometry_1.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = 0.0002839577,
	.box = { BOHR(20.0), BOHR(20.0), BOHR(20.0) },
	.opts = {
		.terms = EFP_TERM_ELEC,
		.elec_damp = EFP_ELEC_DAMP_OFF,
		.enable_pbc = 1,
		.enable_cutoff = 1,
		.swf_cutoff = BOHR(6.0)
	}
};

DEFINE_TEST(elec_1c)
