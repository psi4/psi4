/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/

#ifndef _psi_bin_oeprop_globals_h_
#define _psi_bin_oeprop_globals_h_

#define MAXMP 3         /* Up to octopole moments can be computed */
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define PFACCUTOFF 1.0E-15	/* If AO OPDM matrix element < PFACCUTOFF
				   don't compute the integral */
#define MAXFACT 100     /* Half-size of the double factorial array df*/
#define MFOPDMLOC 34	/* Location of AO density in master file */
#define ADOTB_ORTHOGONAL 1.0E-8 /* Dot product less than this -> vectors are orthogonal */

#define PRINTOPDMLEVEL 3	/* Printing level to print out density matrix */
#define PRINTNMOLEVEL 4		/* Printing level to print out natural MOs */
#define PRINTDIPMLEVEL 3	/* Printing level to print out dipole moment
				   integrals */
#define PRINTOVERLAPLEVEL 4	/* Printing level to print out overlap matrix */
#define PRINTZVECTORLEVEL 5	/* Printing level to print out Z-vector */
#define PRINTBASSETLEVEL 4	/* Printing level to print out basis set
				   information */
#define PRINTTASKPARAMLEVEL 1	/* Printing level to print out options and
				   calculations parameters */
#define PRINTAOPOPLEVEL 2	/* Printing level to print out Mulliken AO
				   population matrix */
#define PRINTCCOEFFLEVEL 5	/* Printing level to print out coupling
				   coefficient vectors */
#define PRINTOCCUPLEVEL 5	/* Printing level to print out occupation
				   vectors */
#define PRINTDIPCOMPLEVEL 2	/* Printing level to print out contributions
				   to electric dipole moment */
#define PRINTDARWINCOMPLEVEL 2	/* Printing level to print out contributions
				   to Darwin term */
#define MAXDENSGRAD 3.0		/* Cutoff on the density gradient */

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

namespace psi { namespace oeprop {

EXTERN int *ioff;
EXTERN double df[MAXFACT*2];

	/* Calculation constants */

EXTERN int nrho; /* number of densities to analyze */
EXTERN int irho; /* index of current density */
EXTERN char **opdm_lbl;
EXTERN char **opdm_a_lbl;
EXTERN char **opdm_b_lbl;
EXTERN int natom, natom3, openmos, openirrs, iopen, nsym, nirreps, charge;
EXTERN int nbfao, natri, nbfso, nstri, nshell, nprim;
EXTERN int nmo;
EXTERN int *clsdpi;		/* Array of numbers of closed shells per irrep */
EXTERN int *openpi;		/* Same for open shells */
EXTERN int *orbspi;		/* (molecular) orbitals per irrep */
EXTERN int *sopi;               /* symmetry orbitals per irrep */
EXTERN int *sprim, *snuc, *stype, *snumg, *sloc, *sloc_new;	/* See documentation for libchkpt library */
EXTERN int lmax;	/* Highest angular momentum of Gaussians in the basis set */
EXTERN char **irr_labs;		/* Irrep labels */
EXTERN std::string title;		/* Calculation title */
EXTERN double **geom;		/* Cartesian geometry */
EXTERN double *zvals, *exps;	/* Nuclear charges and exponents */
EXTERN double *contr;		/* Contraction coefficients */
EXTERN double **scf_evec_so, **scf_evec_ao;	/* SCF eigenvector in two forms */
EXTERN double *scf_evals;	/* SCF eigenvalues */
EXTERN double **usotao;				/* SO to AO transformation matrix */
EXTERN double **usotbf;     /* SO to BF (5d/7f) transformation matrix */


	/* Calculation options */

EXTERN bool read_opdm;		/* Flag for reading density from disk */
EXTERN int opdm_file;		/* Density matrix file number */	
EXTERN std::string opdm_basis;	/* In what basis a onepdm to be read */
EXTERN std::string opdm_format;	/* Format of onepdm file (lower triangle or square) */
EXTERN bool transdens;           /* Read transition densities vs reg dens?  */
EXTERN bool asymm_opdm;		/* Flag for symmetrization of opdm read in from a file */
EXTERN bool wrtnos;		/* Flag for writing NOs to file30 */
EXTERN int print_lvl;		/* Overall printing level */
EXTERN bool print_nos;           /* Print natural orbitals? */
EXTERN bool spin_prop;		/* Write dipole moment to ASCII file ? */
EXTERN int corr;		/* Correlation correction to the first-order 
				   properties flag */
EXTERN int mpmax;		/* Compute up to electric mpmax-tuple moment */
EXTERN double mp_ref_xyz[3];	/* Coordinates of the reference point for elec. mult. moment 
				   calculations */
EXTERN int mp_ref;		/* Code of the reference point for elec. 
				   multipole moment calculations :
				   0 = default (currently - center of mass)
				   1 = center of mass (COM)
				   2 = origin of the space coordinate system
				   3 = center of electronic charge
				   4 = center of nuclear charge
				   5 = center of net charge 
				   If MP_REF_XYZ is specified - MP_REF
				   keyword is set to -1 */
EXTERN double Lm_ref_xyz[3];	/* Coordiantes of the reference point for elec. angular momentum 
				   calculations */
EXTERN int wrt_dipmom;		/* Flag for writing dipole moments into dipmom.dat */
EXTERN bool nuc_esp;		/* Flag for computing electrostatic properties (such as 
				   electrostatic potential, electric field, and 
				   electric field gradient) at the nuclei */
EXTERN int grid;		/* 0 = compute nothing
				   1 = compute electrostatic potential on a 2d grid
				   2 = compute electron density on a 2d grid
				   3 = compute electron density gradient on a 2d grid
				   4 = compute Laplacian of the electron density on a 2d grid
				   5 = evaluate MOs on a 3d grid
				   6 = evaluate density on a 3d grid */
EXTERN int num_mos_to_plot;     /* total number of MOs to plot */
EXTERN int *mos_to_plot;        /* index of the MOs to plot (in Pitzer order) */
EXTERN int grid3d;              /* 1 if this is a 3d grid, 0 otherwise */
EXTERN std::string grid_format;       /* output format for the grid data */
EXTERN double grid_origin[3];   /* Origin of the grid coordinate system as specified by user, then
				   the origin of the grid rectangle/box in the reference system */
EXTERN double grid_unit_x[3];   /* Unit vectors of an intermediate coordinate system in which the grid
				   rectangle/box will be defined (expressed in the reference system) */
EXTERN double grid_unit_y[3];
EXTERN double grid_unit_z[3];   /* Only used if a 3d grid is specified */
EXTERN double grid_xy0[2];	/* A vertex of the grid rectangle (if a 2d grid is used) */
EXTERN double grid_xy1[2];	/* A diagonally opposite vertex of the grid rectangle (if a 2d grid is used) */
EXTERN double grid_xyz0[3];	/* A vertex of the grid box (if a 3d grid is used) */
EXTERN double grid_xyz1[3];	/* A diagonally opposite vertex of the grid box (if a 3d grid is used) */
EXTERN int nix,niy,niz;		/* Number of intervals along x, y, and z (if 3d grid) axes */
EXTERN double grid_step_x[3];   /* Unit grid step along the x axis of the grid rectangle/box (in terms of the reference system) */
EXTERN double grid_step_y[3];
EXTERN double grid_step_z[3];   /* Only used if a 3d grid is specified */
EXTERN double grid_zmin;	/* Lower limit on displayed values of a scalar property on a 2d grid */
EXTERN double grid_zmax;	/* Upper limit on displayed values of a scalar property on a 2d grid */
EXTERN int edgrad_logscale;     /* If non-zero then to use logarithmic scaling of the electron density gradient */
EXTERN int zvec_file;		/* Z-vector file number */
EXTERN int delete_zvec;		/* Whether to delete the Z-vector file */


	/* Density matrices */

EXTERN double *Ptot, *Pspin;


	/* Integral intermediates */

EXTERN double ***MIX, ***MIY, ***MIZ;   /* "Boxes" of (lmax+1) by (lmax+1) by 3 
    (up to octopole moment is supported) of elementary moment integrals
    defined as in Obara and Saika paper JCP 84 (1986) 3963 (eqs. A3-A8). 
    Center C coincides with P. */
EXTERN double ***AI0, ***AIX, ***AIY, ***AIZ;    /* Elementary integrals of A-operator (see Obara and Saika paper) */
EXTERN double ***AIXX, ***AIYY, ***AIZZ, ***AIXY, ***AIXZ, ***AIYZ;


	/* Arrays of exponents of x, y, and z in basis functions and "relative"
	   normalization constants (see initialize.c) */

EXTERN int **xpow_bf, **ypow_bf, **zpow_bf;
EXTERN double **norm_bf;


	/* Integral matrices */

EXTERN double *S;


	/* Z-vector in AO basis : 
           Zmn = Sum[Cim*Cjn*Zij,{i,0,nbfso-1},{j,0,i}] 
           where Zij is Z-vector in MO basis, 
	   Cim - SCF eigenvector in AO basis. */

EXTERN double **zvec;


	/* Results of the calculation */

EXTERN double *qnet;		/* Net atomic charge */
EXTERN double dx,dy,dz;	/* Components of dipole moment */
EXTERN double dx_e,dy_e,dz_e,dx_n,dy_n,dz_n;
EXTERN double dxcc,dycc,dzcc;	/* Correlation correction to the dipole moment */
EXTERN double dtot;		/* Total dipole moment */
EXTERN double qxx,qxy,qxz,qyy,qyz,qzz;	/* Components of the quadrupole moment */
EXTERN double qxxcc,qyycc,qzzcc,
	      qxycc,qxzcc,qyzcc;	/* Correlation corrections to the
					   quadrupole moment */
EXTERN double *qvals;	/* Principal values of the quadrupole moment tensor */
EXTERN double **qvecs;	/* Principal axis of the quadr. tensor */
EXTERN double exp_x2, exp_y2, exp_z2;	/* Expectation values of x^2, y^2, z^2
					   operators */
EXTERN double *MOXX, *MOYY, *MOZZ;	/* Orbital spatial extents */
EXTERN double oxxx,oyyy,ozzz,oxxy,oxxz,	/* Components of the octopole moment */
              oxyy,oyyz,oxzz,oyzz,oxyz;
EXTERN double oxxxcc,oyyycc,ozzzcc,oxxycc,oxxzcc,	/* Corrections */
              oxyycc,oyyzcc,oxzzcc,oyzzcc,oxyzcc;
EXTERN double Lx,Ly,Lz;	        /* Components of electronic angular momentum L */
EXTERN double Lx2,Ly2,Lz2;      /* Components of square electronic angular momentum L^2 */
EXTERN double *phi;		/* Electrostatic potential at the nuclei */
EXTERN double *ex, *ey, *ez;	/* Components of electric field at the nuclei */
EXTERN double *dexx, *deyy, *dezz, *dexy, *dexz, *deyz;	/* Gradients of E */
EXTERN double *ahfsxx, *ahfsyy, *ahfszz,
              *ahfsxy, *ahfsxz, *ahfsyz;  /* Dipole-dipole contributions to the hyperfine coupling constants */
EXTERN double *edens, *sdens;   /* Electron and spin densities at the nuclei */
EXTERN double massveloc, darw;  /* First-order relativistic corrections to the energy */
EXTERN double *darw_per_atom; /* atomic contributions to one-electron Darwin term */
EXTERN double **grid_pts;	/* Scalar property (electrostatic potential, electron density, density gradient magnitude, Laplacian)
				   values on a 2d grid */
EXTERN double **grid_vecX, **grid_vecY, **grid_vecZ;    /* Components of vector properties (density gradient,
							      principal values of the Hessian) on the grid */
EXTERN double ****grid3d_pts;	/* Values of properties on a 3d grid:
				 -first dimension runs over each property (MOs in MO case, density,
				  density gradient components, etc.)
				 -second, third, and fourth dimensions run over x, y, and z indices of grid points */
EXTERN double **nmo_so;		/* Natural orbitals in the SO basis */
EXTERN double **nmo_ao;		/* Natural orbitals in the AO basis */

EXTERN int **connectivity;    /* A matrix of flags whether a bong exists between two atoms */
EXTERN std::string wfn;             /* wavefunction type */
EXTERN std::string ref;             /* reference type */
EXTERN int update_energy_with_MVD; /* update energy in file 1? for SCF_MVD optimizations */
EXTERN double fine_structure_alpha; /* multiply relativisitic terms by this */
EXTERN bool QED_darwin; /* scale atomic contributions to get QED correction */

}} // namespace psi::oeprop

#endif // header guard
