#include "test_common.h"

static const char *files[] = {
	"./fraglib/h2o.efp",
	"./fraglib/c6h6.efp",
	"./fraglib/nh3.efp",
	 NULL
};

static const char *names[] = {
	"H2O_L",
	"C6H6_L",
	"NH3_L",
	 NULL
};

static const double frag_xyzabc[] = {
	BOHR(-1.6), BOHR( 4.7), BOHR( 1.4), -1.3,  0.1,  7.0,
	BOHR( 0.4), BOHR(-0.9), BOHR(-0.7),  2.3,  1.6, -2.3,
	BOHR(-3.5), BOHR(-2.0), BOHR(-0.7),  0.0,  2.2,  2.7
};

static const double qm_znuc[] = {
	1.0, 8.0, 2.0, 1.0
};

static const double qm_xyz[] = {
	BOHR( 3.2), BOHR( 1.8), BOHR(-2.3),
	BOHR(-2.9), BOHR(-6.2), BOHR(-2.5),
	BOHR( 5.0), BOHR( 4.3), BOHR( 0.2),
	BOHR( 4.9), BOHR( 0.0), BOHR( 4.7)
};

static enum efp_result get_electron_density_field(int n_pt,
		UNUSED const double *xyz, double *field,
		UNUSED void *user_data)
{
	/* no electrons */
	memset(field, 0, n_pt * 3 * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = frag_xyzabc,
	.n_ptc = ARRAY_SIZE(qm_znuc),
	.ptc_charges = qm_znuc,
	.ptc_xyz = qm_xyz,
	.ref_energy = -0.0787829370,
	.electron_density_field_fn = get_electron_density_field,
	.electron_density_field_user_data = NULL,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR |
			 EFP_TERM_AI_ELEC | EFP_TERM_AI_POL,
		.elec_damp = EFP_ELEC_DAMP_OVERLAP,
		.disp_damp = EFP_DISP_DAMP_TT
	}
};

DEFINE_TEST(qm_1)
