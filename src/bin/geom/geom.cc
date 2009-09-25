/*! \defgroup GEOM geom: Compute and print geometrical parameters */

/*! 
** \file
** \brief Compute geometrical parameters
** \ingroup GEOM
**
** Program calculate geometrical parameters from data provided
** in format of PSI's FILE11.DAT or GEOM.DAT
**
** Reads the last section in file11.dat, or the first section
** in geom.dat, as it should.  Should also read in coordinates
** that are in a format similar to geom.dat but not quite.
**
** Does interatomic distances, bond angles, out-of-plane (oop)
** angles, torsional angles.
**
** Also, center of mass, principal moments of intertia, 
** rotational constants, if the data has been read in file11.dat
** format (which contains atomic numbers).
**
** Created by C. David Sherrill in April 1993
**
** Modified by CDS:
**  7/93 to read geometry from geom.dat 
**  9/94 for more flexible reading of geom.dat and for -h flag support
** 11/94 to increase MAXATOM to 60 (might need to do Buckeyball!)
**       Also completely remove dependence on all libraries except 
**       libciomr (& thus also libipv1).  Included function 
**       fill_sym_matrix() locally to remove dependence on my math library.
**       Changed function malloc_check() to malloc_ck() due to conflict with
**       a libciomr function of the same name.
**  9/95 to print atomic coordinates in angstroms also
** Modified by MLL:
**  7/95 to do isotopomers and use new masses header file.
** Modified by CDS:
**  5/98 to be able to read input in angstroms
**  7/99 to read Q-Chem format
** 
*/

/* include's */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <masses.h>
#include <physconst.h>


/* define's */
#define MAX_LINE 132
#define MAXATOM 100
#define AN_TOL 0.0001
#define ITOL 0.100
#define MIN_I 0.00001
#define PRINT_DIST 4.00  /* max pair distance for printing (bohr) */

namespace psi { namespace geom {

/* global variables */
static char line1[MAX_LINE];
char line2[MAX_LINE];
int print_all_dist = 0; /* flag for printing all parameters, however far
                         * apart individual pairs of atoms are */
double print_dist = PRINT_DIST;
extern "C" {
FILE *outfile;                 /* output file pointer */
FILE *infile; 
char *psi_file_prefix;
}

/* declare functions in this file */
void malloc_ck(void *array, const char *mesg);
void prf_abort(FILE *file, const char *mesg);
int read_file11(char* label, int* natom, double* energy, double** X, double** Y, double** Z, double** AN, FILE* fpo);
void print_file11(char* label, int natom, double energy, double* X, double* Y, double* Z, double* AN, FILE* fpo) ;
int read_geom(int maxlines, int* natom, double** X, double** Y, double** Z, char* fname);
void print_geom(int natom, double *X, double *Y, double *Z, FILE *fpo);
int read_aces_geom(char *fname, int *natom, double **X, double **Y,
                   double **Z, double **AN, FILE *fpo);
void print_aces_geom(int natom, double *X, double *Y, double *Z, double *AN,
                     const char **labels, FILE *fpo);
int read_qchem_geom(char *fname, int *natom, double **X, double **Y,
                    double **Z, double **AN, FILE *fpo);
int read_xyz_geom(char *fname, int *natom, double **X, double **Y,
                  double **Z, double **AN, FILE *fpo);
void calc_distances(int natom, double *X, double *Y, double *Z, double ***R);
void calc_unit_vectors(int natom, double *X, double *Y, double *Z,
                       double **R, double E[MAXATOM][MAXATOM][3]);
void calc_bond_angles(int natom, double E[MAXATOM][MAXATOM][3],
                      double BondAngles[MAXATOM][MAXATOM][MAXATOM],
                      double **R, FILE *fpo);
void calc_oop_angles(int natom, double E[MAXATOM][MAXATOM][3],
                     double BondAngles[MAXATOM][MAXATOM][MAXATOM],
                     double **R, FILE *fpo);
void calc_tors_angles(int natom, double E[MAXATOM][MAXATOM][3],
                      double BondAngles[MAXATOM][MAXATOM][MAXATOM],
                      double **R, FILE *fpo);
double dot_vect(double *a, double *b, int n);
void cross_vect(double *a, double *b, double *c);
void center_of_mass(int natom, double *M, double *X, double *Y, double *Z,
                    double CM[3]);
void calc_relative_coords(int natom, double *X, double *Y, double *Z,
                          double **XR, double **YR, double **ZR, double CM[3]);
void calc_itensor(int natom, double *M, double *X, double *Y, double *Z,
                  double **I);
void calc_molec_type(int natom, double Ia, double Ib, double Ic, FILE *fpo);
void calc_rot_constants(double Ia, double Ib, double Ic, FILE *fpo);
void calc_mass_analysis(int natom, double *M, double *X, double *Y, 
			double *Z, FILE *fpo);
void fill_sym_matrix(double **A, int size);
}} // namespace psi::geom

int main(int argc, char* argv[])
{
  using namespace psi::geom;
  int natom;                     /* number of atoms */
  double energy;                 /* energy from file11 */
  double *X, *Y, *Z;             /* pointers to arrays of Cartesian coords */
  double **Tmat;                 /* temp matrix of Cartesians in angstroms */
  double *M;                     /* pointer to mass array */
  double *AN;                    /* pointer to atomic number array */
  double **R;                    /* interatomic distance matrix */
  double **RA;                   /* R, in angstroms */
  double E[MAXATOM][MAXATOM][3]; /* unit vector matrix */
  char label[MAX_LINE];          /* label on file11 */
  char outfname[] = "geom.out";  /* file to output to */
  int i, j, k;                   /* loop variables */
  double BondAngles[MAXATOM][MAXATOM][MAXATOM];  /* 3d arr for bond angles */
  int do_oop = 0;                /* default not to do out of plane angles */
  int ang_in = 0;                /* input is in angstroms */
  char geom_file[50];            /* filename for geometry file */
  int errcod = 0;                /* necessary for input parsing */
  int num_array;                 /* number of arrays of isotopes */
  char tmpstr[100];              /* tmp string for geom isotopes array */
  int tmpi = -1;                 /* tmp stor for casting atomic nums to int */
  int num_unparsed=0;            /* number of unparsed cmd line args */
  char *argv_unparsed[100];      /* pointers to unparsed cmd line args */
  double tval;                   
  enum GeomFormat { 
    PSI_FILE11,
    PSI_GEOM,
    ACES2,
    QCHEM,
    XYZ
  } InputType;


  /* Initialize and zero data */ 
  for (i=0; i<MAXATOM; i++) { 
    for (j=0; j<MAXATOM; j++) {
      for (k=0; k<MAXATOM; k++) {
	if ( k < 3) E[i][j][k] = 0.0;
	BondAngles[i][j][k] = 0.0;
      }
    }
  }
    
  InputType = PSI_FILE11;
   
  /* get command line arguments */
  strcpy(geom_file, "geom.dat");
  for (i=1, num_unparsed=0; i<argc; i++) {
    if (strcmp("-g", argv[i]) == 0) {
      InputType = PSI_GEOM;
      if ( (i+1 < argc) && (strchr(argv[i+1], '-') == NULL) ) {
	strcpy(geom_file, argv[i+1]) ;
        i++;
      }
    }
    else if (strcmp("-aces", argv[i]) == 0) {
      InputType = ACES2;
      if ( (i+1 < argc) && (strchr(argv[i+1], '-') == NULL) ) {
	strcpy(geom_file, argv[i+1]) ;
        i++;
      }
    }
    else if (strcmp("-qchem", argv[i]) == 0) {
      InputType = QCHEM;
      ang_in = 1;  /* QCHEM format is in Angstrom */
      if ( (i+1 < argc) && (strchr(argv[i+1], '-') == NULL) ) {
	strcpy(geom_file, argv[i+1]) ;
        i++;
      }
    }
    else if (strcmp("-xyz", argv[i]) == 0) {
      InputType = XYZ;
      ang_in = 1;  /* XYZ format is in Angstrom */
      if ( (i+1 < argc) && (strchr(argv[i+1], '-') == NULL) ) {
	strcpy(geom_file, argv[i+1]);
        i++;
      }
    }
    else if (strcmp("-oop", argv[i]) == 0) do_oop = 1;
    else if (strcmp("-a", argv[i]) == 0) print_all_dist = 1;
    else if (strcmp("-angstroms", argv[i]) == 0) ang_in = 1;
    else if (strcmp("-angstrom", argv[i]) == 0) ang_in = 1;
    else if (strcmp("-d", argv[i]) == 0) { 
      if (sscanf(argv[i+1], "%lf", &print_dist)!=1) print_dist = PRINT_DIST;
      i++;
    }
    else if (strcmp("-h", argv[i]) == 0) {
      fprintf(stdout,
	      "geom: options available\n");
      fprintf(stdout,
	      " -h         = help (this list)\n");
      fprintf(stdout,
	      " -g [fname] = read geometry in Cartesians from geom.dat or\n");
      fprintf(stdout,
	      "              specified file.\n");
      fprintf(stdout,
	      " -aces [fname]  = read geom in Aces output format\n");
      fprintf(stdout,
	      " -qchem [fname] = read geom in QCHEM output format\n");
      fprintf(stdout,
	      " -xyz [fname] = read geom in XYZ output format\n");
      fprintf(stdout,
	      "                 from geom.dat or specified file.\n");
      fprintf(stdout,
	      " -oop       = print out-of-plane angles\n");
      fprintf(stdout,
	      " -d dist    = print params only for pairs < dist apart\n");
      fprintf(stdout,
	      "              (default distance is 4.0 bohr).\n");
      fprintf(stdout,
	      " -a         = print params for ALL pair distances\n");
      fprintf(stdout,
	      " -angstrom  = input is in angstroms\n");
      exit(0);
    }
    else {
      argv_unparsed[num_unparsed++] = argv[i];
    }
  }

  /* open files */
  errcod = psi_start(&infile,&outfile,&psi_file_prefix,num_unparsed,argv_unparsed,0);
  ip_cwk_add(":GEOM");
  /* done inside psi_start now
  ffile(&outfile,outfname,0);
  ffile(&infile,"input.dat",2); 
  */
  tstart(outfile);

  /* check input.dat for running parameters */
  errcod = ip_boolean("READ_GEOM",&i,0);
  if ((errcod == IPE_OK) && i==1) InputType = PSI_GEOM;
  errcod = ip_boolean("ACES",&i,0);
  if ((errcod == IPE_OK) && i==1) InputType = ACES2;
  errcod = ip_boolean("QCHEM",&i,0);
  if ((errcod == IPE_OK) && i==1) InputType = QCHEM;
  errcod = ip_boolean("XYZ",&i,0);
  if ((errcod == IPE_OK) && i==1) InputType = XYZ;
  errcod = ip_boolean("PRINT_ALL_DIST",&print_all_dist,0);
  errcod = ip_boolean("DO_OOP",&do_oop,0);
  errcod = ip_boolean("ANGSTROM",&ang_in,0);
  errcod = ip_boolean("ANGSTROMS",&ang_in,0);
  errcod = ip_data("PRINT_DIST","%lf",&print_dist,0);

  /* print program identification */
  fprintf(outfile, "          ***********************************\n");
  fprintf(outfile, "          **   GEOMETRY ANALYSIS PROGRAM   **\n");
  fprintf(outfile, "          **     David Sherrill,  1993     **\n");   
  fprintf(outfile, "          ***********************************\n");
  fprintf(outfile, "\n");
  fprintf(outfile, "   Use the -h flag to list available options\n\n");


  /* first, read in the data */
  switch (InputType) {
  case PSI_FILE11:
    if (!read_file11(label, &natom, &energy, &X, &Y, &Z, &AN, outfile)) {
      prf_abort(outfile, "geom: trouble reading file11\n");
    }
    print_file11(label, natom, energy, X, Y, Z, AN, outfile);
    break;
  case PSI_GEOM:
    if (!read_geom(MAXATOM, &natom, &X, &Y, &Z, geom_file)) {
      sprintf(line1, "geom: trouble reading geom file %s\n", geom_file);
      prf_abort(outfile, line1);
    }
    fprintf(outfile, "GEOMETRY DATA INPUT\n") ;
    fprintf(outfile, "Number of atoms = %d\n", natom) ;
    fprintf(outfile, "Cartesian coordinates (%s):\n", ang_in ? 
	    "angstroms" : "bohr");
    print_geom(natom, X, Y, Z, outfile);
    break;
  case ACES2:
    if (!read_aces_geom(geom_file, &natom, &X, &Y, &Z, &AN, outfile)) {
      sprintf(line1, "geom: trouble reading geom file %s\n", geom_file);
      prf_abort(outfile, line1);
    }
    fprintf(outfile, "DATA FROM INPUT\n"); 
    fprintf(outfile, "Number of atoms = %d\n", natom); 
    fprintf(outfile, "\nCartesian coordinates (bohr) :\n"); 
    print_aces_geom(natom, X, Y, Z, AN, atomic_labels, outfile);
    break;
  case QCHEM:
    if (!read_qchem_geom(geom_file, &natom, &X, &Y, &Z, &AN, outfile)) {
      sprintf(line1, "geom: trouble reading geom file %s\n", geom_file);
      prf_abort(outfile, line1);
    }
    fprintf(outfile, "GEOMETRY DATA INPUT\n") ;
    fprintf(outfile, "Number of atoms = %d\n", natom) ;
    fprintf(outfile, "Cartesian coordinates (%s):\n", ang_in ? 
	    "angstroms" : "bohr");
    print_geom(natom, X, Y, Z, outfile);
    break;
  case XYZ:
    if (!read_xyz_geom(geom_file, &natom, &X, &Y, &Z, &AN, outfile)) {
      sprintf(line1, "geom: trouble reading geom file %s\n", geom_file);
      prf_abort(outfile, line1);
    }
    fprintf(outfile, "GEOMETRY DATA INPUT\n") ;
    fprintf(outfile, "Number of atoms = %d\n", natom) ;
    fprintf(outfile, "Cartesian coordinates (%s):\n", ang_in ? 
	    "angstroms" : "bohr");
    print_aces_geom(natom, X, Y, Z, AN, atomic_labels, outfile);
    break;
  default:
    prf_abort(outfile, "Couldn't figure out your input type!\n");
  }
  
  tval = _bohr2angstroms;
  if (ang_in) tval = 1.0 / tval;
  
  Tmat = init_matrix(3,natom);
  for (i=0; i<natom; i++) {
    Tmat[0][i] = X[i] * tval;
    Tmat[1][i] = Y[i] * tval;
    Tmat[2][i] = Z[i] * tval;
  }
  
  fprintf(outfile, "Cartesian coordinates (%s):\n", ang_in ?
	  "bohr" : "angstroms");
  if (InputType == PSI_GEOM || InputType == QCHEM)
    print_geom(natom, Tmat[0], Tmat[1], Tmat[2], outfile);
  else 
    print_aces_geom(natom, Tmat[0], Tmat[1], Tmat[2], AN, atomic_labels, 
                    outfile);
  
  if (ang_in) {
    for (i=0; i<natom; i++) {
      X[i] = Tmat[0][i];
      Y[i] = Tmat[1][i];
      Z[i] = Tmat[2][i];
    }
  }
  free_matrix(Tmat, 3);
  
  if (natom > MAXATOM) {
    sprintf(line1, "geom: maximum number of atoms %d exceeded\n", MAXATOM) ;
    prf_abort(outfile, line1) ;
  }
  
  /* calculate all interatomic distances */
  calc_distances(natom, X, Y, Z, &R) ;
  fprintf(outfile, "\nInteratomic distance matrix (bohr)\n") ;
  print_mat(R, natom, natom, outfile) ;
  fprintf(outfile, "\n") ;
  
  /* also print interatomic distances in angstroms */
  RA = init_matrix(natom, natom) ;
  for (i=0; i<natom; i++) {
    for (j=0; j<natom; j++) {
      RA[i][j] = R[i][j] * _bohr2angstroms ;
    }
  }
  fprintf(outfile, "\nInteratomic distance matrix (angstroms)\n") ;
  print_mat(RA, natom, natom, outfile) ;
  fprintf(outfile, 
	  "\nGeometrical Parameters for pairs less than %.2f bohr apart\n", 
	  print_dist) ;
  fprintf(outfile, "(for all parameters use -a flag)\n") ;
  fflush(outfile) ;
  
  
  
  /* calculate all bond angles */
  calc_unit_vectors(natom, X, Y, Z, R, E) ;
  calc_bond_angles(natom, E, BondAngles, R, outfile) ;
  fflush(outfile) ;
  
  /* calculate out of plane angles */
  if ( (natom > 3) && do_oop) {
    calc_oop_angles(natom, E, BondAngles, R, outfile) ;
  }
  fflush(outfile) ;
  
  /* calculate the torsional angles */
  if (natom > 3) {
    calc_tors_angles(natom, E, BondAngles, R, outfile) ;
    fflush(outfile) ;
  }
  
  /* get masses, center-of-mass, rotational constants */
  
  M = init_array(natom); 
  
  errcod = ip_count("ISOTOPES", &num_array, 0);
  if (errcod != IPE_OK) {
    if (InputType != PSI_GEOM && InputType != QCHEM) {
      for (i=0 ; i<natom ; i++) {
	tmpi = (int)AN[i]; 
	M[i] = an2masses[tmpi]; 
      }
      calc_mass_analysis(natom, M, X, Y, Z, outfile);
    }
    
    /* close files */
    fprintf(outfile, "\n");
    tstop(outfile);
    psi_stop(infile,outfile,psi_file_prefix);
    exit(0);
  }
  
  
  for (i=0 ; i<num_array ; i++) {
    if (num_array > 1) 
      fprintf(outfile, "\n****** Analysis for isotopomer %d ******\n\n", i);
    errcod = ip_count("ISOTOPES", &k, 1, i) ; 
    if (errcod != IPE_OK || k != natom) {
      printf("GEOM: trouble parsing isotopes array %d\n", i) ; 
      fprintf(outfile, "\n") ;
      tstop(outfile);
      psi_stop(infile,outfile,psi_file_prefix);
      exit(0) ;
    } 
    
    for (j=0 ; j<natom ; j++) {
      errcod = ip_data("ISOTOPES","%s", tmpstr,2,i,j) ;
      if (isdigit(tmpstr[0])) {
	sscanf(tmpstr, "%lf", &M[j]) ;
      }
      else {
	for (k=0 ; k<LAST_MASS_INDEX; k++) {
	  if (strcmp(tmpstr, mass_labels[k]) == 0) {
	    M[j] = atomic_masses[k] ;
	    break ;
	  }
	}
	if (k == LAST_MASS_INDEX) {
	  printf("GEOM: Problem with element %d of array %d", j+1, i+1);
	  printf("in ISOTOPES section of input.dat\n");
	  exit(0) ;
	}
      }
    }
    calc_mass_analysis(natom, M, X, Y, Z, outfile);
  } /* end parsing of ISOTOPES keyword */
  
  /* close files */
  fprintf(outfile, "\n");
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}


namespace psi { namespace geom {

/*
** CALC_MASS_ANALYSIS(): This function reads in an array of the
**    mass of each atom, calculates the center of mass, translates
**    the coordinates of each atom relative to the center of mass,
**    forms the moment of inertia tensor and diagonalizes it,
**    sorts the principal moments of inertia in order Ia<Ib<Ic,
**    and determines the molecular type and rotational constants.
*/
void calc_mass_analysis(int natom, double *M, double *X, double *Y, 
			double *Z, FILE *fpo)
{
  int i, j ;
  double **I ;                   /* moment of inertia tensor */
  double **I_evecs ;             /* eigenvectors of I tensor */
  double *I_evals ;              /* eigenvalues of I tensor */
  double Ia, Ib, Ic ;            /* principal moments of inertia */
  double *XR, *YR, *ZR ;         /* coords relative to C of M */
  double CM[3] ;                 /* center of mass coordinates */
  
  
  /* allocate the memory */
  I = init_matrix(3, 3) ;
  I_evecs = init_matrix(3, 3) ;
  I_evals = (double *) malloc (3 * sizeof(double)) ;
  malloc_ck((void *)I_evals, "geom: trouble allocating I eigenval array\n");
  
  
  /* get the center of mass */
  center_of_mass(natom, M, X, Y, Z, CM) ;
  fprintf(fpo,"\nCenter of Mass Coords (bohr):\n   %12.7lf %12.7lf %12.7lf\n",
	  CM[0], CM[1], CM[2]);
  fprintf(fpo,
        "\nCenter of Mass Coords (angstroms):\n   %12.7lf %12.7lf %12.7lf\n",
	  CM[0]*_bohr2angstroms, CM[1]*_bohr2angstroms, CM[2]*_bohr2angstroms);
  fflush(fpo);
  
  /* translate the atomic coords rel to center of mass */
  calc_relative_coords(natom, X, Y, Z, &XR, &YR, &ZR, CM) ;
  fprintf(fpo, "\nCoordinates relative to center of mass (bohr):\n") ;
  for (i=0; i<natom; i++) {
    fprintf(fpo, "     %10.6lf  %15.10lf  %15.10lf  %15.10lf\n", 
	    M[i], XR[i], YR[i], ZR[i]) ;
  }
  
  /* form the moment of inertia tensor I */
  calc_itensor(natom, M, XR, YR, ZR, I) ;
  fprintf(fpo, "\nMoment of Inertia Tensor") ;
  print_mat(I, 3, 3, fpo) ;
  fflush(fpo) ;
  
  /* diagonalize the moment of inertia tensor */
  sq_rsp(3, 3, I, I_evals, 1, I_evecs, 1.0E-14);
  fprintf(fpo, "\nInertia tensor eigenvectors\n") ;
  eivout(I_evecs, I_evals, 3, 3, fpo) ;
  
  /* sort out which is Ia, Ib, Ic */
  Ia = I_evals[0] ; j = 0 ;
  for (i=0; i<3; i++) 
    if (I_evals[i] < Ia) { Ia = I_evals[i] ; j = i ; }
  I_evals[j] = 9999999.9 ; /* make sure it can't be chosen again */
  Ib = I_evals[0] ; j = 0 ;
  for (i=0; i<3; i++) 
    if (I_evals[i] < Ib) { Ib = I_evals[i] ; j = i ; }
  I_evals[j] = 9999999.9 ; /* make sure it can't be chosen again */
  Ic = I_evals[0] ;
  for (i=0; i<3; i++) 
    if (I_evals[i] < Ic) Ic = I_evals[i] ; 
  fprintf(fpo, "\nPrincipal Moments of Inertia (a.u.):\n") ;
  fprintf(fpo, "   Ia = %12.6f\n", Ia) ;
  fprintf(fpo, "   Ib = %12.6f\n", Ib) ;
  fprintf(fpo, "   Ic = %12.6f\n", Ic) ;
  fprintf(fpo, "\nPrincipal Moments of Inertia \n") ;
  fprintf(fpo, "            amu A**2       g cm**2\n") ;
  fprintf(fpo, "   Ia   %12.6f   %12.6g\n",
	  Ia * _bohr2angstroms * _bohr2angstroms, 
	  Ia * _bohr2cm * _bohr2cm * _amu2g) ;
  fprintf(fpo, "   Ib   %12.6f   %12.6g\n",
	  Ib * _bohr2angstroms * _bohr2angstroms, 
	  Ib * _bohr2cm * _bohr2cm * _amu2g) ;
  fprintf(fpo, "   Ic   %12.6f   %12.6g\n", 
       Ic * _bohr2angstroms * _bohr2angstroms,
	  Ic * _bohr2cm * _bohr2cm * _amu2g) ;
  
  /* now convert everything to amu A**2 */
  Ia = Ia * _bohr2angstroms * _bohr2angstroms ;
  Ib = Ib * _bohr2angstroms * _bohr2angstroms ;
  Ic = Ic * _bohr2angstroms * _bohr2angstroms ;
  
  /* calculate the molecular type (linear, spherical top, etc.) */
  calc_molec_type(natom, Ia, Ib, Ic, fpo) ;
  
  /* calculate the rotational constants */
  calc_rot_constants(Ia, Ib, Ic, fpo) ;
  fflush (fpo) ;
}



/*
** CALC_DISTANCES(): This function, given the number of atoms, and the
**   Cartesian coordinates of the atoms, calculates the matrix of
**   interatomic distances and stores it in R (which it also allocates)
*/
void calc_distances(int natom, double *X, double *Y, double *Z, double ***R) 
{
  int i, j ;
  double delx, dely, delz, dist;
  
  /* Now calculate all possible interatomic distances R(ij) */
  *R = init_matrix(natom, natom);
  /* more efficient if we calculate only lower triangle of R(ij) */
  for (i=0; i<natom; i++) {
    for (j=0; j<=i; j++) {
      delx = X[i] - X[j];
      dely = Y[i] - Y[j];
      delz = Z[i] - Z[j];
      dist = delx * delx + dely * dely + delz * delz;
      dist = sqrt(dist);
      (*R)[i][j] = dist;
    }
  }
  fill_sym_matrix(*R, natom);
  
}


/*
** CALC_UNIT_VECTORS(): Function calculates unit vectors between each
**    pair of atoms and stores them in E.
*/
void calc_unit_vectors(int natom, double *X, double *Y, double *Z, 
		       double **R, double E[MAXATOM][MAXATOM][3])  
{
  int i, j ;
  
  for (i=0; i<natom; i++) {
    for (j=0; j<natom; j++) {
      if (i != j) {
	E[i][j][0] = -(X[i] - X[j]) / R[i][j] ;
	E[i][j][1] = -(Y[i] - Y[j]) / R[i][j] ;
	E[i][j][2] = -(Z[i] - Z[j]) / R[i][j] ;
      }
      else {
	E[i][j][0] = 0.0 ;
	E[i][j][1] = 0.0 ;
	E[i][j][2] = 0.0 ;
      }
    }
  }
}


/*
** CALC_BOND_ANGLES(): Function calculates all non-redundant bond angles
**   between all combinations of three atoms, and writes to file fpo.
**   Note: the first index of BondAngles is always smaller; i.e.
**   angle 5-3-1 is stored only as 1-3-5
*/
void calc_bond_angles(int natom, double E[MAXATOM][MAXATOM][3], 
		      double BondAngles[MAXATOM][MAXATOM][MAXATOM], 
		      double **R, FILE *fpo)
{
  int i, j, k;
  double dotprod;
  double angle;
  
   fprintf(fpo, "\nBond Angles:\n");
   for (i=0; i<natom; i++) {
     for (j=0; j<natom; j++) {
       for (k=i; k<natom; k++) {
	 if ( (i != j) && (i != k) && (j != k) ) {
	   dotprod = dot_vect(E[j][i], E[j][k], 3) ;
	   if (dotprod > 1.00000) angle = 0.0000 ;
	   else if (dotprod < -1.00000) angle = _pi ;
	   else angle = acos(dotprod) ;
	   BondAngles[i][j][k] = angle ;
	   if (((R[i][j] < print_dist) && (R[j][k] < print_dist)) ||
	       print_all_dist) {
	     fprintf(fpo, "%2d-%2d-%2d    %14.8lf\n", i+1, j+1, k+1, 
                     (angle*180.00/(double)_pi)) ;
	   }
	 }
       }
     }
   }
}

/*
** MALLOC_CK(): Function checks to see if pointer supposedly pointing
**    to a recently malloc'ed memory space actually points to anything.
**    If not, assume error is fatal, print message, and abort
**
** Arguments: 
**       array    =  pointer to memory supposedly allocated (casted to void)
**       mesg     =  character string containing error message 
*/
void malloc_ck(void *array, const char *mesg)
{
  if (array == NULL) {
    fprintf(stderr, "%s", mesg) ;
    exit(0) ;
  }
}


/*
** PRF_ABORT(): Function prints error message to file, closes that file,
**     and aborts.
**
** Arguments:
**       file    =  file pointer to write error to 
**       mesg    =  character string containing error message 
*/
void prf_abort(FILE *file, const char *mesg)
{
  printf("%s", mesg);
  fprintf(file, "%s", mesg);
  fclose(file);
  exit(0);
}


/*
** DOT_VECT(): Function takes dot product of two linear arrrays (a and b)
**    of size n.
**
** Returns: the dot product (double)
*/
double dot_vect(double *a, double *b, int n)
{
  register int i;
  double tval = 0.0;
  
  for (i=0; i<n; i++) 
    tval += a[i] * b[i];
  
  return(tval);     
}

/*
   Function to compute the cross product of 2 Cartesian vectors
*/

void cross_prod(double *v1, double *v2, double *out)
{
   out[0] =     v1[1]*v2[2]-v1[2]*v2[1];
   out[1] =    -v1[0]*v2[2]+v1[2]*v2[0];
   out[2] =     v1[0]*v2[1]-v1[1]*v2[0];
   return;
}

/*
** CALC_TORS_ANGLES(): This function calculates the torsional angles
**   and writes output to a file.  Define 'impossible' angles (e.g.
**   torsional angles where 3 atoms are collinear) as 0.
** 
** Arguments:
**       natom      = number of atoms
**       BondAngles = 3D array, gives angle between atoms i,j,k
**       fpo        = file pointer for output file
**       E          = 2D array of unit vectors. 
**                    E[i][j] is unit vector between atoms i and j
**       R          = Interatomic distance (bohr) matrix
*/
void calc_tors_angles(int natom, double E[MAXATOM][MAXATOM][3], 
		      double BondAngles[MAXATOM][MAXATOM][MAXATOM], 
		      double **R, FILE *fpo)
{
  int i, j, k, l, xyz ;
  double cross1[3], cross2[3], cross3[3] ;
  double tval, angle, phi2, phi3 ;
  double sign, norm3;
  
  fprintf(fpo, "\nTorsional Angles:\n") ;
  for (i=0; i<natom; i++) {
    for (j=0; j<natom; j++) {
      for (k=0; k<natom; k++) {
	for (l=i; l<natom; l++) {
	  if ( (i!=j) && (i!=k) & (i!=l) && (j!=k) && (j!=l) && (k!=l)) {
	    cross_vect(E[i][j], E[j][k], cross1) ;
	    cross_vect(E[j][k], E[k][l], cross2) ;
	    if (i < k) phi2 = BondAngles[i][j][k] ;
	    else phi2 = BondAngles[k][j][i] ;
	    if (j < l) phi3 = BondAngles[j][k][l] ;
	    else phi3 = BondAngles[l][k][j] ;
	    tval = dot_vect(cross1, cross2, 3) ;
	    if ((sin(phi2) > 0.00001) && (sin(phi3) > 0.00001)) {
	      tval /= sin(phi2) ; 
	      tval /= sin(phi3) ;
	    }
	    else tval = 2.0 ;
	    if (tval > 0.99999) angle = 0.0000 ;
	    else if (tval < -0.99999) angle = _pi ;
	    else angle = acos(tval) ;

	    /* compute the sign */
	    cross_prod(cross1, cross2, cross3);
	    norm3 = sqrt(dot_vect(cross3,cross3,3));
	    sign = 1.0;
	    if (fabs(norm3) > 0.00001) {
	      for(xyz=0; xyz<3; ++xyz)
		cross3[xyz] *= 1.0/norm3;
	      tval = dot_vect(cross3, E[j][k],3);
	      if (tval < 0.0)
		sign = -1.0;
	    }

	    if ( ((R[i][j] < print_dist) && (R[j][k] < print_dist) &&
		  (R[k][l] < print_dist)) || print_all_dist) {
	      fprintf(fpo, "%2d-%2d-%2d-%2d    %13.8lf\n",
		      i+1, j+1, k+1, l+1, sign * angle * 180.0/_pi) ;
	    }
	  }
	}
      }
    }
  }
}


/*
** CALC_OOP_ANGLES(): This function calculates the out of plane angles
**   and writes output to a file.
** 
** Arguments:
**       natom      = number of atoms
**       BondAngles = 3D array, gives angle between atoms i,j,k
**       fpo        = file pointer for output file
**       E          = 2D array of unit vectors. 
**                    E[i][j] is unit vector between atoms i and j
**       R          = interatomic distance (bohr) matrix
*/
void calc_oop_angles(int natom, double E[MAXATOM][MAXATOM][3], 
		     double BondAngles[MAXATOM][MAXATOM][MAXATOM], 
		     double **R, FILE *fpo)
{
  int i, j, k, l;
  double tval1, tval2, angle;
  double crossprod[3];
  
  fprintf(fpo, "\nOut of Plane Angles:\n") ;
  fprintf(fpo, "  Out of plane angle is the angle formed by the vector 1-4\n");
  fprintf(fpo, "  and the plane defined by atoms 2, 3, and 4\n");
  
  for (i=0; i<natom; i++) {
    for (j=0; j<natom; j++) {
      for (k=j; k<natom; k++) {
	for (l=0; l<natom; l++) {
	  if ( (i!=j) && (i!=k) & (i!=l) && (j!=k) && (j!=l) && (k!=l)) {
	    cross_vect(E[l][j], E[l][k], crossprod) ;
	    tval1 = dot_vect(crossprod, E[l][i], 3) ;
	    /* make sure j<k when BondAngles is accessed */
	    if (j < k)  tval2 = BondAngles[j][l][k] ; 
	    else tval2 = BondAngles[k][l][j] ; 
	    if (sin(tval2) > 0.0001) tval1 = tval1 / sin(tval2) ;
	    else tval1 = 0.0 ;
	    if (tval1 > 1.000000) angle = (_pi / 2.00000) ;
	    else if (tval1 < -1.000000) angle = -(_pi / 2.0000) ;
	    else angle = asin(tval1) ;
	    if ( ( (R[i][j] < print_dist) && (R[j][k] < print_dist) &&
		   (R[k][l] < print_dist) ) || print_all_dist ) {
	      fprintf(fpo, "%2d-%2d-%2d-%2d    %13.8lf\n", 
		      i+1, j+1, k+1, l+1, angle * 180.0/(double)_pi) ;
	    }
	  }
	}
      }
    }
  }
}


/*
** CROSS_VECT(): Take cross product of two vectors (3-element arrays only)
**    
** Arguments: a, b = two vectors to cross (each 3 element array of doubles)
**             c   = vector to hold results of a x b
*/
void cross_vect(double *a, double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1] ;
  c[1] = a[2] * b[0] - a[0] * b[2] ;
  c[2] = a[0] * b[1] - a[1] * b[0] ;
}


/*
** CENTER_OF_MASS(): Find the center of mass of the molecule, given
**   the mass and cartesian coordinate arrays
**
** Arguments: 
**     natom   =  number of atoms
**         M   =  mass array
**     X,Y,Z   =  cartesian coordinates arrays
**        CM   =  3-element array (double) for center of mass coords
*/
void center_of_mass(int natom, double *M, double *X, double *Y, double *Z, 
		    double CM[3])
{
  double tot_mass = 0.0;
  double sumx = 0.0, sumy = 0.0, sumz = 0.0;
  int i;
  
  for (i=0; i<natom; i++) {
    tot_mass += M[i];
    sumx += M[i] * X[i];
    sumy += M[i] * Y[i];
    sumz += M[i] * Z[i];
  }
  CM[0] = sumx / tot_mass;
  CM[1] = sumy / tot_mass;
  CM[2] = sumz / tot_mass;
  
}


/*  The calc_mass_array routine is commented out. 

** CALC_MASS_ARRAY(): This function takes the array holding the atomic
**    number of each atom and finds the atomic mass from a datafile
**    masses.dat
**
** Returns: 1 if success, 0 if fail

int calc_mass_array(int natom, double *AN, double **M, FILE *fpo)
{
   FILE *fpi ;
   int i ;
   double an, am ;

   *M = (double *) malloc (natom * sizeof(double)) ;
   malloc_ck((void *) *M, "(calc_mass_array): can't alloc mass array\n") ;
   if ((fpi = fopen(MASSFILE, "r")) == NULL) {
      printf("geom: error opening mass file %s\n", MASSFILE);
      exit(0);
      }
   fprintf(fpo, "\nNuclear masses (amu): \n") ;
   for (i=0; i<natom; i++) {
      rewind(fpi) ;
      while (fscanf(fpi, "%lf %lf\n", &an, &am) == 2) {
         if (fabs(an - AN[i]) < AN_TOL) {
            (*M)[i] = am ;
            fprintf(fpo, "%2d  %f  %.7lf\n", i+1, AN[i], am) ;
            break ;
            }
         }
      if (feof(fpi) || ferror(fpi)) {
         prf_abort(fpo, "(calc_mass_array): Can't find atomic number data\n");
         }
      }
   fprintf(fpo, "\n") ;
   return(1) ;
}
  
*/


/*
** CALC_RELATIVE_COORDS(): Form arrays of coordinates relative to
**    center of mass.
*/
void calc_relative_coords(int natom, double *X, double *Y, double *Z, 
			  double **XR, double **YR, double **ZR, double CM[3]) 
{
  int i;

  *XR = (double *) malloc (natom * sizeof(double)) ;
  malloc_ck((void *) *XR, 
	    "(calc_relative_coords): Trouble allocating array\n");
  *YR = (double *) malloc (natom * sizeof(double)) ;
  malloc_ck((void *) *YR, 
	    "(calc_relative_coords): Trouble allocating array\n");
  *ZR = (double *) malloc (natom * sizeof(double)) ;
  malloc_ck((void *) *ZR, 
	    "(calc_relative_coords): Trouble allocating array\n");
  
  for (i=0; i<natom; i++) {
    (*XR)[i] = X[i] - CM[0];
    (*YR)[i] = Y[i] - CM[1];
    (*ZR)[i] = Z[i] - CM[2];
  }
}


/*
** CALC_ITENSOR(): Form the moment of inertia tensor, given the mass
**    array and cartesian arrays.
**
** Arguments: 
**     natom   =  number of atoms
**         M   =  mass array
**     X,Y,Z   =  cartesian arrays
**         I   =  pointer to moment of inertia tensor (allocated here)
*/
void calc_itensor(int natom, double *M, double *X, double *Y, double *Z, 
		  double **I) 
{
  int i, j;
  
  /* zero elements */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      I[i][j] = 0.0;
    }
  }
  
  /* do diagonal ones first */
  for (i=0; i<natom; i++) {
    I[0][0] += M[i] * (Y[i] * Y[i] + Z[i] * Z[i]) ;
    I[1][1] += M[i] * (X[i] * X[i] + Z[i] * Z[i]) ;
    I[2][2] += M[i] * (X[i] * X[i] + Y[i] * Y[i]) ;
  }
  
  /* do off-diagonal ones */
  for (i=0; i<natom; i++) {
    I[1][0] -= M[i] * (Y[i] * X[i]) ; /* I(y,x) */
    I[2][0] -= M[i] * (Z[i] * X[i]) ; /* I(z,x) */
    I[2][1] -= M[i] * (Z[i] * Y[i]) ; /* I(z,y) */
  }
  I[0][1] = I[1][0] ;     /* I(x,y) */
  I[0][2] = I[2][0] ;     /* I(x,z) */
  I[1][2] = I[2][1] ;     /* I(y,z) */
  
}


/*
** CALC_MOLEC_TYPE(): Function determines if molecule is in one of several
**   common classes (linear, diatomic, spherical top, etc.) 
**
** Arguments: 
**        natom  = number of atoms
**        Ia, Ib, Ic = Principal moments of inertia
**        fpo    = file pointer to write results to
*/
void calc_molec_type(int natom, double Ia, double Ib, double Ic, FILE *fpo)
{
  double diff_ab, diff_bc;
  
  fprintf(fpo, "\nMolecular type: ");
  
  if (natom == 2) {
    fprintf(fpo, "diatomic\n");
    return;
  }
  diff_ab = Ia - Ib;  if (diff_ab < 0.0) diff_ab = -diff_ab;
  diff_bc = Ib - Ic;  if (diff_bc < 0.0) diff_bc = -diff_bc;
  
  if ( (Ia < ITOL) && (diff_bc < ITOL) ) {
    fprintf(fpo, "linear\n");
    return;
  }
  else if ( (diff_ab < ITOL) && (diff_bc < ITOL) ) {
    fprintf(fpo, "spherical top\n");
    return;
  }
  else if (diff_ab < ITOL) {
    fprintf(fpo, "oblate top\n");
    return;
  }
  else if (diff_bc < ITOL) {
    fprintf(fpo, "prolate top\n");
    return;
  }
  else fprintf(fpo, "asymmetric top\n");
  
}

      
/* 
** CALC_ROT_CONSTANTS(): Calculates the rotational constants
**   of a molecule, given the principal moments of inertia
**   Ia, Ib, and Ic.  Writes output to file fpo.
**
** n.b. Function expects moments of inertia in amu * angstroms^2 
*/
void calc_rot_constants(double Ia, double Ib, double Ic, FILE *fpo) 
{
  double A1, A2, B1, B2, C1, C2, tval1, tval2 ;
  int noA = 0 ;
  
  /* first convert amu * angstroms^2 to kg*m^2 */
  if (Ia < MIN_I) noA = 1 ;
  Ia = Ia * _amu2kg ; Ia *= (double) 1E-20 ;
  Ib = Ib * _amu2kg ; Ib *= (double) 1E-20 ;
  Ic = Ic * _amu2kg ; Ic *= (double) 1E-20 ;
  
  /* now do wavenumbers */
  tval1 = _h / (8.0 * _pi * _pi * _c * 100.0) ;
  if (!noA) A1 = tval1 / Ia ; 
  B1 = tval1 / Ib; C1 = tval1 / Ic ;
  
  /* do MHz */
  tval2 = _h / (8000000.0 * _pi * _pi) ;
  if (!noA) A2 = tval2 / Ia ; 
  B2 = tval2 / Ib; C2 = tval2 / Ic ;
  
  fprintf(fpo, "\nRotational constants: \n") ;
  
  /* print out A, B, C: different formats for different units */
  fprintf(fpo, "             cm**-1         MHz\n") ;
  if (noA) { /* if no A constant, just print B and C */
    fprintf(fpo, "     A = %10.4f       %10.2f\n", B1, B2) ;
    fprintf(fpo, "     B = %10.4f       %10.2f\n", C1, C2) ;
  }
  else {
    fprintf(fpo, "     A = %10.4f       %10.2f\n", A1, A2) ;
    fprintf(fpo, "     B = %10.4f       %10.2f\n", B1, B2) ;
    fprintf(fpo, "     C = %10.4f       %10.2f\n", C1, C2) ;
  }
}

void fill_sym_matrix(double **A, int size)
{
  double **row, *col;
  int rc, cc;
  
  row = A;
  for (rc = 0; rc < (size-1); rc++) {
    col = *row;
    for (cc = 0; cc < size; cc++) {
      if (cc > rc) {
	*col = A[cc][rc];
      }
      col++;
    }
    row++;
  }
}


/*
** label2an
**
** This function returns the atomic number corresponding to a given mass
** label, or -1 if not found.
**
** C. David Sherrill
** July 1999
*/
int label2an(char *label)
{
  int i, j, k;
  char p, q;

  for (i=0; i<LAST_ATOMIC_INDEX; i++) {
    k = strlen(label);
    for (j=0; j < k; j++) {
       p = label[j];
       q = atomic_labels[i][j];
       if (toupper(p) != toupper(q)) break; 
    }
    if (j == k) return(i);
  }

  /* couldn't find the label */
  return(-1);

}

}} // namespace psi::geom
