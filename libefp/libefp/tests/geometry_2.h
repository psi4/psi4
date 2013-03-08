static const char *files[] = {
	"./fraglib/h2o.efp",
	"./fraglib/nh3.efp",
	 NULL
};

static const char *names[] = {
	"H2O_L",
	"NH3_L",
	"H2O_L",
	"H2O_L",
	"NH3_L",
	 NULL
};

static const double xyzabc[] = { /* some random geometry */
	BOHR(-1.0), BOHR( 3.7), BOHR( 0.4), -1.3,  0.0,  7.0,
	BOHR( 0.4), BOHR(-0.9), BOHR(-0.7),  4.0,  1.6, -2.3,
	BOHR( 1.7), BOHR( 2.0), BOHR( 3.3), -1.2, -2.0,  6.2,
	BOHR( 0.0), BOHR( 3.9), BOHR(-3.4),  1.3,  5.2, -3.0,
	BOHR(-3.5), BOHR( 0.0), BOHR(-0.7),  0.0, -2.7,  2.7
};
