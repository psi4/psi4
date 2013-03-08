static const char *files[] = {
	"./fraglib/acetone.efp",
	"./fraglib/c2h5oh.efp",
	"./fraglib/c6h6.efp",
	"./fraglib/ccl4.efp",
	"./fraglib/ch3oh.efp",
	"./fraglib/ch4.efp",
	"./fraglib/cl2.efp",
	"./fraglib/dcm.efp",
	"./fraglib/dmso.efp",
	"./fraglib/h2.efp",
	"./fraglib/h2o.efp",
	"./fraglib/nh3.efp",
	 NULL
};

static const char *names[] = {
	"ACETONE_L",
	"C2H5OH_L",
	"C6H6_L",
	"CCL4_L",
	"CH3OH_L",
	"CH4_L",
	"CL2_L",
	"DCM_L",
	"DMSO_L",
	"H2_L",
	"H2O_L",
	"NH3_L",
	 NULL
};

static const double xyzabc[] = { /* some random geometry */
	BOHR( 0.0), BOHR( 0.0), BOHR(0.0), 0.0, 0.2, 0.3,
	BOHR( 7.0), BOHR( 0.0), BOHR(0.0), 0.0, 2.0, 3.7,
	BOHR(14.0), BOHR( 0.0), BOHR(0.0), 3.1, 0.8, 2.0,
	BOHR(21.0), BOHR( 0.0), BOHR(0.0), 0.0, 8.0, 0.0,
	BOHR( 0.0), BOHR( 6.0), BOHR(0.0), 0.7, 2.0, 1.0,
	BOHR( 7.0), BOHR( 6.0), BOHR(0.0), 0.6, 0.0, 4.7,
	BOHR(14.0), BOHR( 6.0), BOHR(0.0), 0.0, 0.0, 0.3,
	BOHR(21.0), BOHR( 6.0), BOHR(0.0), 0.0, 0.4, 0.3,
	BOHR( 0.0), BOHR(12.0), BOHR(0.0), 0.8, 0.0, 0.0,
	BOHR( 7.0), BOHR(12.0), BOHR(0.0), 8.0, 0.7, 0.8,
	BOHR(14.0), BOHR(12.0), BOHR(0.0), 0.0, 0.0, 0.0,
	BOHR(21.0), BOHR(12.0), BOHR(0.0), 0.0, 2.0, 0.0
};
