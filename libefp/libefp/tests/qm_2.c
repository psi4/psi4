#include "test_common.h"

static const char *files[] = {
	"./fraglib/ch3oh.efp",
	"./fraglib/dmso.efp",
	"./fraglib/dcm.efp",
	"./fraglib/acetone.efp",
	 NULL
};

static const char *names[] = {
	"CH3OH_L",
	"DMSO_L",
	"DMSO_L",
	"ACETONE_L",
	"DCM_L",
	"ACETONE_L",
	"ACETONE_L",
	 NULL
};

static const double frag_xyzabc[] = {
	BOHR( 0.0), BOHR(-1.0), BOHR( 0.0),  0.0,  1.1,  2.0,
	BOHR(-5.0), BOHR(12.0), BOHR(-0.0),  3.0,  0.2,  5.0,
	BOHR(-0.0), BOHR(-3.0), BOHR( 5.0),  6.0,  2.3,  8.0,
	BOHR(-5.0), BOHR(-4.0), BOHR(-5.0),  9.0,  0.4,  1.0,
	BOHR(-9.0), BOHR(-5.0), BOHR( 1.0),  2.0,  1.5,  4.0,
	BOHR( 7.0), BOHR(-2.0), BOHR(11.0),  5.0,  0.6,  7.0,
	BOHR(-9.0), BOHR(-7.0), BOHR(-9.0),  8.0,  2.7,  0.0
};

static const double qm_znuc[] = {
	7.0, 10.0, 1.0, 2.0, 2.0
};

static const double qm_xyz[] = {
	BOHR( 4.0), BOHR( 5.0), BOHR( 5.0),
	BOHR( 8.0), BOHR( 5.0), BOHR( 6.0),
	BOHR( 5.0), BOHR( 8.0), BOHR( 5.0),
	BOHR( 8.0), BOHR( 9.0), BOHR( 8.0),
	BOHR( 5.0), BOHR( 8.0), BOHR( 8.0)
};

static enum efp_result get_electron_density_field(int n_pt,
		const double *xyz, double *field, void *user_data)
{
	(void)xyz;
	(void)user_data;

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
	.ref_energy = -0.2314262632,
	.electron_density_field_fn = get_electron_density_field,
	.electron_density_field_user_data = NULL,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR |
			 EFP_TERM_AI_ELEC | EFP_TERM_AI_POL,
		.elec_damp = EFP_ELEC_DAMP_SCREEN,
		.disp_damp = EFP_DISP_DAMP_OVERLAP
	}
};

DEFINE_TEST(qm_2)
