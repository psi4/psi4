#include "test_common.h"
#include "geometry_6.h"

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = xyzabc,
	.ref_energy = -0.0051253344,
	.box = { BOHR(15.0), BOHR(15.0), BOHR(15.0) },
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR,
		.disp_damp = EFP_DISP_DAMP_TT,
		.elec_damp = EFP_ELEC_DAMP_OVERLAP,
		.pol_damp = EFP_POL_DAMP_TT,
		.enable_pbc = 1,
		.enable_cutoff = 1,
		.swf_cutoff = BOHR(5.0)
	}
};

DEFINE_TEST(total_6b)
